%% DEMO: Truss Optimization

clear; clc; close all

%generate a random trimesh in a box
V = [0 0; 4 0; 4 1; 0 1];
n = 40;
[V]=evenlySampleCurve(V,n,'linear',1);
pointSpacing = 4/n;                 %node spacing in mesh
stdP=0.5*pointSpacing*ones(1,2);    %randomness in mesh
[Ft,Vt,boundaryNodes]=regionTriMeshRand2D({V},pointSpacing,stdP,0,0);
DT = triangulation(Ft,Vt);
LI = edges(DT);

%find the index of the fixed nodes
fixed = find(Vt(:,1)==0);

%collect some loaded nodes at the center of the top boundary
forced = find(sum((Vt-[4 0.5]).^2,2)<0.05);

%define load
load = 0.01*[0 -1]';

%solve for displacements
[D,C] = LDS_Bar_Solver(DT,1,fixed,forced,load);

%compute stresses carried by each member
S = deformedStressPlot(DT,D,1,1);

plotV(Vt(fixed,:),'b.','markersize',20)


