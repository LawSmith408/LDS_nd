%% DEMO: Truss Optimization

clear; clc; close all

plotme = true;

%generate a random trimesh in a box
V = [-1 0; 1 0; 1 1; -1 1];
n = 50;
[V]=evenlySampleCurve(V,n,'linear',1);
[Ft,Vt,boundaryNodes]=regionTriMeshRand2D({V},4/n,[1/n 1/n],0,0);
DT = triangulation(Ft,Vt);
LI = edges(DT);

%find the index of the fixed nodes
lowerRight = find(sum((Vt-[ 1 0]).^2,2)<0.05);
lowerLeft  = find(sum((Vt-[-1 0]).^2,2)<0.05);
fixed = [lowerRight; lowerLeft];

%collect some loaded nodes at the center of the top boundary
forced = find(sum((Vt-[0 1]).^2,2)<0.05);

%define load
load = 0.01*[0 -1]';

%initialize materials
rho = 0.5*ones(size(edges(DT),1),1);

Emin = 1e-3;
Emax = 1e3;
rhomin = 1e-1;
rhomax = 1-1e-1;
p=3;

for i = 1:8

%compute per-element elasticity
eMat = Emin + rho.^p*(Emax-Emin);

%solve for displacements
[D,C] = LDS_Bar_Solver(DT,eMat,fixed,forced,load);

%compute stresses carried by each member
S = deformedStressPlot(DT,D,eMat,0);

cla
rhoPlot(DT,rho); 
plotV(Vt(fixed,:),'b.','markersize',20)
plotV(Vt(forced,:),'r.','markersize',20)
drawnow();

%compute per-element compliance sensitivities
%reference: https://link.springer.com/article/10.1007/s001580050176
for j = 1:length(eMat)
    V = Vt(LI(j,:),:);
    [k,L] = barLocalStiffness(V);
    u = D(LI(j,:),:)';
    dcdx(j) = p*rho(j)^(p-1)*(Emax-Emin)*u(:)'*eMat(j)*k*u(:);
    dvdx(j) = L; %total lattice volume is Sigma (L_i * rho_i)
end
lambda = -1;
nu = 0.5;
Be = -dcdx./(lambda*dvdx);
move = rho.*Be'.^nu;
rho = rho + move;

% %update density
% rho = rho.*(rescale(dcdx(:)).^0.05);

%saturate
rho(rho<rhomin) = rhomin;
rho(rho>rhomax) = rhomax;

end

% cla
% deformedDensityPlot(DT,D,rescale(rho)); 
% plotV(Vt(fixed,:),'b.','markersize',20)
% plotV(Vt(forced,:),'r.','markersize',20)
% drawnow();

function [k,l] = barLocalStiffness(V)

%V is a 2xdim matrix of the [X Y (Z)] coords of a pair of points

DV = diff(V);           %subtract the endpoints
l = sqrt(sum(DV.^2));   %compute the length of this member
C = DV/l;               %compute the cosine angles Cx Cy Cz
L = C'*C;               %stiffness submatrix lambda
k = [L -L; -L L];       %full stiffness matrix

end
