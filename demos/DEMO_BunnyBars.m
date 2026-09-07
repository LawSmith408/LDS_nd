clear; clc; close all
setupLDS    %put src/, demos/ and input_meshes/ on the MATLAB path

%% load geometry
% load small_bunny_mesh2D.mat
load small_bunny_mesh2D_vfine.mat
%load small_bunny_mesh3D.mat
%load small_bunny_mesh3D_fine.mat


%define boundary conditions
NC = DT.Points;
fixed = find(NC(:,1)<-40);
forced = find(NC(:,1)>15);

%solve system
[D,C] = LDS_Bar_Solver(DT,1,fixed,forced);

%% Generate Output Plot
cFigure
patch('faces',DT.ConnectivityList,'vertices',NC+D,'FaceColor','w','handlevisibility','off');hold on
plotV(NC(fixed,:)+D(fixed,:),'r.','markersize',20)
plotV(NC(forced,:)+D(forced,:),'b.','markersize',20)
legend('fixed','loaded'); axis equal
if size(NC,2)>2; axisGeom; end; gdrawnow()