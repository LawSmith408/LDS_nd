%% DEMO: SIMP topology optimization of a bar lattice
% Lawrence Smith | lasm4254@colorado.edu
%
% Minimize the compliance of a truss subject to a volume constraint, using
% SIMP material interpolation and an optimality-criteria (OC) update.
%
%   design variable   rho(e)  relative cross-sectional area of bar e, in [rhomin,1]
%   interpolation     A(rho) = Emin + rho^p*(E0-Emin)          (SIMP)
%   objective         c(rho) = F'*u = sum_e (A_e/L_e)*dL_e^2   (compliance)
%   constraint        sum_e rho_e*L_e <= volfrac*sum_e L_e     (material budget)
%
% The whole optimization is ~30 lines because LDS_Bar_Solver already accepts
% a per-element stiffness vector, and the compliance sensitivity of a bar
% has a closed form that needs no element loop:
%
%   dc/drho_e = -p*rho_e^(p-1)*(E0-Emin)/L_e * dL_e^2
%
% where dL_e is the axial elongation of bar e. That falls straight out of
% u_e'*k0_e*u_e = dL_e^2 for the unit bar stiffness matrix k0.
%
% Reference: Bendsoe & Sigmund, "Material interpolation schemes in topology
% optimization", https://link.springer.com/article/10.1007/s001580050176

clear; clc; close all
setupLDS    %put src/, demos/ and input_meshes/ on the MATLAB path

%% Problem setup: simply supported span, point load at the top centre
V = [-1 0; 1 0; 1 1; -1 1];
n = 50;
[V] = evenlySampleCurve(V,n,'linear',1);
[Ft,Vt] = regionTriMeshRand2D({V},4/n,[1/n 1/n],0,0);
DT = triangulation(Ft,Vt);
LI = edges(DT);                             %every edge is a bar

%pin the two lower corners
fixed = [find(sum((Vt-[ 1 0]).^2,2)<0.05); find(sum((Vt-[-1 0]).^2,2)<0.05)];

%load the top centre downward
forced  = find(sum((Vt-[0 1]).^2,2)<0.05);
loadVec = 0.01*[0 -1]';

%% Precompute element geometry (fixed throughout: only rho changes)
dV   = Vt(LI(:,2),:)-Vt(LI(:,1),:);         %edge vectors
L0   = vecnorm(dV,2,2);                     %edge lengths
Cdir = dV./L0;                              %unit direction cosines

%% Optimization parameters
E0      = 1;        %solid stiffness
Emin    = 1e-6;     %void stiffness (keeps K non-singular)
p       = 3;        %SIMP penalization exponent
volfrac = 0.30;     %fraction of the ground structure's material we may keep
rhomin  = 1e-2;     %lower bound on density
move    = 0.2;      %OC move limit
eta     = 0.5;      %OC damping exponent
maxIter = 60;
tol     = 5e-3;     %stop when max density change falls below this

volTarget = volfrac*sum(L0);
rho = volfrac*ones(size(LI,1),1);           %start on the volume constraint

%The system is deliberately ill conditioned (Emin:E0 spans 1e6), so let the
%solver quietly fall back to its direct solve instead of warning each pass.
warning('off','LDS:pcgFailed');
cleanup = onCleanup(@() warning('on','LDS:pcgFailed'));

%% Optimization loop
history = nan(maxIter,2);   %[compliance, volume fraction]
change  = inf;
iter    = 0;

figure('Position',[385 200 1044 344]);

while change > tol && iter < maxIter
    iter = iter+1;

    %SIMP interpolation. LDS_Bar_Solver scales the *unit* bar stiffness by
    %eMat, so eMat plays the role of EA/L -- divide by length here.
    A    = Emin + rho.^p*(E0-Emin);
    eMat = A./L0;

    %solve the equilibrium problem
    [D,c] = LDS_Bar_Solver(DT,eMat,fixed,forced,loadVec);

    %axial elongation of every bar, vectorized
    dL = sum(Cdir.*(D(LI(:,2),:)-D(LI(:,1),:)),2);

    %compliance and volume sensitivities
    dc = -p*rho.^(p-1)*(E0-Emin)./L0 .* dL.^2;      %always negative
    dv = L0;                                        %d/drho of sum(rho.*L)

    %optimality-criteria update, with lambda bisected onto the volume constraint
    rhoNew = ocUpdate(rho,dc,dv,volTarget,move,rhomin,eta);
    change = max(abs(rhoNew-rho));
    rho    = rhoNew;

    %sanity check on the sensitivity derivation: the compliance returned by
    %the solver must equal the sum of element strain energies computed from
    %the elongations. If these disagree, dc below is wrong.
    if iter==1
        fprintf('check: solver c = %.6e, sum(eMat.*dL.^2) = %.6e\n',...
            c,sum(eMat.*dL.^2));
    end

    history(iter,:) = [c sum(rho.*L0)/sum(L0)];
    fprintf('it %3i | c %.4e | vol %.3f | change %.4f\n',iter,c,history(iter,2),change);

    %% Plot progress
    subplot(1,2,1); cla
    rhoPlot(DT,rho); set(gca,'CLim',[0 1])   %CLim, not clim(), for pre-R2022a
    plotV(Vt(fixed,:), 'b.','markersize',10)
    plotV(Vt(forced,:),'r.','markersize',10)
    title(sprintf('iteration %i',iter)); axis off

    subplot(1,2,2); cla
    yyaxis left;  plot(1:iter,history(1:iter,1),'-','linewidth',1.5); ylabel('compliance')
    yyaxis right; plot(1:iter,history(1:iter,2),'-','linewidth',1.5); ylabel('volume fraction')
    ylim([0 1]); xlabel('iteration'); xlim([1 max(2,iter)]); grid on
    title('convergence')

    drawnow();
end

fprintf('Converged in %i iterations. Compliance %.4e (started %.4e).\n',...
    iter,history(iter,1),history(1,1));

%% Optimality-criteria update
function rhoNew = ocUpdate(rho,dc,dv,volTarget,move,rhomin,eta)
%Bisect the Lagrange multiplier lambda until the updated design sits exactly
%on the volume constraint. The updated density is monotonically decreasing
%in lambda, which is what makes the bisection valid.

l1 = 0; l2 = 1e9;
while (l2-l1)/(l1+l2) > 1e-9 && l2 > 1e-30    %second test guards l1==0
    lmid   = 0.5*(l1+l2);
    B      = -dc./(lmid*dv);                        %positive: dc<0, dv>0
    rhoNew = max(rhomin, max(rho-move, min(1, min(rho+move, rho.*B.^eta))));
    if sum(rhoNew.*dv) > volTarget
        l1 = lmid;      %too much material, raise the price
    else
        l2 = lmid;
    end
end

end
