%% DEMO: Truss Optimization
%Author - Lawrence Smith

clear; clc; close all

%generate a random trimesh in a box
V = [-1 0; 1 0; 1 1; -1 1];
n =100;
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
% 
% patch('faces',Ft,'vertices',Vt,'facecolor','w'); hold on
% plotV(Vt(forced,:),'b.','markersize',30)
% plotV(Vt(fixed,:),'r.','markersize',30)


%define load
load = 0.01*[0 -1]';

%initialize materials
rho = 0.5*ones(size(edges(DT),1),1);

set(gcf,'Position',[385 200 1044 344])
% rhoPlot(DT,rho); 

Emin = 1e-3;    Emax = 1;
rhomin = 1e-1; rhomax = 1-1e-1;
p=2;

%try just running 10 iterations
for i = 1:10

%compute per-element elasticity
eMat = Emin + rho.^p*(Emax-Emin);

%solve for displacements
[D,C] = LDS_Bar_Solver(DT,eMat,fixed,forced,load);

%compute stresses carried by each member
S = deformedStressPlot(DT,D,eMat,0);

%plot element densities
cla
subplot(1,2,1)
rhoPlot(DT,rho); 
%show boundary conditions
plotV(Vt(fixed,:),'b.','markersize',10)
plotV(Vt(forced,:),'r.','markersize',10)

%plot histogram of element densities
subplot(1,2,2)
histogram(rho,0:0.05:1,'Normalization','probability')
ylim([0 1]);
xlabel('element densities')
ylabel('norm. frequency')

drawnow();

%find which members have a lot of compressive stress
stressedMembers = rescale(S(:))<0.5;

%add some density to the most stressed members
rho(stressedMembers) = rho(stressedMembers)+0.1*rand(nnz(stressedMembers),1);

%subtract some density from the least stressed members
rho(~stressedMembers) = rho(~stressedMembers)-0.1*rand(nnz(~stressedMembers),1);

rho(rho<rhomin) = rhomin;
rho(rho>rhomax) = rhomax;

end

