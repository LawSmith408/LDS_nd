%% Morphing Wing Analysis using LDS Beam Solver
%Lawrence Smith | lasm4254@colorado.edu

%depenencies: GIBBON https://www.gibboncode.org/Installation/

clear; clc; close all

%% Define Load Conditions
%nose cone material properties
E_basematerial = 120e9; %Pa     : modulus of titanium
d_strut = 0.01;         %m      : lattice strut diameter
L_strut = 0.03;         %m      : lattice strut length (average)
Pressure = 0.1;         %atm    : pressure load at tip of cone
nu      = 0.33;         %       : Poisson Ratio

%% Build the Nose Cone geometry 
%nose cone dimensions
r_nose = 0.1;   %m
h_nose = 0.5;  %m

%Define the wing geometry using a parabola
coneFunc = @(x) -(h_nose/r_nose^2)*x.^2 + h_nose;
x = linspace(-1,1,31);
V0 = [[x 0]*r_nose; coneFunc(x*r_nose) 0]'; V = V0;

%generate a triangle mesh inside the region
pointSpacing = 0.02;              %node spacing in mesh
stdP=0.25*pointSpacing*ones(1,2);    %randomness in mesh
[Ft,Vt,boundaryNodes]=regionTriMesh2D({V},pointSpacing,0,0);

%identify nodes which will be constrained
fixed = find(Vt(boundaryNodes(:,1),2)<2*eps);

%identify nodes which will be loaded
forced = find(Vt(boundaryNodes(:,1),2)>0);

%plot the nosecone
patch('Faces',Ft,'Vertices',Vt,'facecolor','none','linewidth',...
    1.5,'displayname','Lattice Struts')
axis equal; hold on
plotV(Vt(fixed,:),'b.','markersize',20,'displayname','Fixed Nodes')
plotV(Vt(forced,:),'r.','markersize',20,'displayname','Loaded Nodes')
set(gcf,'position',[360 175 973 774])
set(gca,'FontName','georgia','fontsize',24)
legend('location','southwest')

%% define the pressure distribution
%determine normal vectors for the cone boundary edge midpoints. We'll do
%this in a discrete sense to stay compatible with non-analytical boundaries
dV0 = V0(1:end-1,:)-V0(2:end,:);  %these vectors 
Vn = vecnormalize(dV0);           %normals of these vectors
LV0 = vecnorm(dV0,2,2);           %lengths of these vectors
Vn = [Vn(:,2) -Vn(:,1)];          %tangents to these vectors

%now average adjacent normal vectors to compute point normals
Vn = (Vn(1:end-1,:)+Vn(2:end,:))*0.5;

%determine the global node IDs for the boundary nodes
for i = 1:length(Vt)
    D = sqrt(sum((Vt(i,:)-V0).^2,2));
    [~,mi] = min(D); map(i) = mi;
end

%define falloff function for scaling pressure load
falloffFunction = @(y) (y/h_nose).^0.5;
Vn = Vn.*falloffFunction(V0(2:end-1,2));
Pload = Pressure*101325*-Vn(map(forced)-1,:)';      %convert to Pa, resolve
PLoad = Pload.*repmat(LV0(2:end-1)',2,1);           %account for area

%now plot the pressure load
scale = 0.1;
Xqv = Vt(forced,1)+Vn(map(forced)-1,1)*scale;
Yqv = Vt(forced,2)+Vn(map(forced)-1,2)*scale;
Uqv = -Vn(map(forced)-1,1)*scale;
Vqv = -Vn(map(forced)-1,2)*scale;
quiver(Xqv,Yqv,Uqv,Vqv,0,'r','linewidth',1.5,'showarrowhead','on',...
    'displayname','Pressure Load')
set(gca,'XColor','none','YColor','none')

%% Analyze the Airfoil using BEAM Elements (support axial and bending loads)
DT = triangulation(Ft,Vt);
LI = edges(DT);
m = size(LI,1); %number of edges
n = size(Vt,1); %number of nodes
Coord = [Vt zeros(size(Vt,1),1)]; %coordinates of nodes
Con = [LI ones(size(LI,1),2)];    %connectivity of each node and whether rotational dofs are held
Re = zeros(n,6);                  %degrees of freedom (start with all free)
Re(fixed,:) = 1;
Load = zeros(n,6); Load(forced,[1 2]) = Pload'; %loading matrixLoad
w = zeros(m,3);
E = E_basematerial*ones(1,m);
G = E./(2*(1+nu));
A = pi*d_strut^2/4*ones(1,m);
Iz = pi*d_strut^4/64*ones(1,m);
Iy = pi*d_strut^4/64*ones(1,m);
J = pi*d_strut^4/32*ones(1,m);
St = zeros(n,6);
be = zeros(1,m);

D=struct('m',m,'n',n,'Coord',Coord','Con',Con','Re',Re',...
    'Load',Load','w',w','E',E','G',G','A',A','Iz',Iz','Iy',Iy',...
    'J',J','St',St','be',be');

[Q,V,R]=MSA(D);

U = V(1:2,:)';

figure
patch('Faces',Ft,'Vertices',Vt,'facecolor','none','edgecolor',[0.8 0.8 0.8],...
    'linewidth',1.5,'displayname','Lattice Struts','handlevisibility','off')
axis equal; hold on
patch('Faces',Ft,'Vertices',Vt+U*10,'facecolor','none','edgecolor','interp',...
    'CData',vecnorm(U*1e3,2,2),'linewidth',2,'displayname','Lattice Struts')
axis equal; hold on
set(gcf,'position',[360 175 973 774])
set(gca,'FontName','georgia','fontsize',16)
set(gca,'XColor','none','YColor','none')
colormap(flipud(brewermap(20,'RdYlBu')))
c = colorbar('location','southoutside');
c.Title.String = 'Displacement [mm]';
title('Uniform Stiffness Cone (disp. magnified 10x)','FontSize',18)

%% Now locally soften a region of the nosecone and re-simulate

%define soft region
c_soft = [r_nose*0.7 h_nose*0.7];

%define a scaling function for how to soften lattice
max_soft_factor = 0.1;  %fraction of full stiffness at center
r_soft= 0.08;      %m
soften = @(r) max_soft_factor+(1-max_soft_factor)./(1+exp(-(r-r_soft)/(0.1*r_soft)));
% r_dummy = linspace(0,4*r_soft);
% plot(r_dummy,soften(r_dummy),'k','linewidth',2)
% xlabel('distance from soft center')
% ylabel('Normalized Stiffness')

%compute the centroid of each lattice element
c_LI = 0.5*(Vt(LI(:,1),:)+Vt(LI(:,2),:));

%compute distance from each centroid to the soft region centerpoint
D_soft = sqrt(sum((c_LI-c_soft).^2,2));
D_soft2= sqrt(sum((Vt-c_soft).^2,2));       %for plotting only

%locally soften selected lattice elements
Evec = E_basematerial*ones(length(LI),1);
Evec_soft = Evec.*soften(D_soft);
Vt_Soft = soften(D_soft2);

figure
patch('Faces',Ft,'Vertices',Vt,'facecolor','none','edgecolor','interp',...
    'CData',Vt_Soft,'linewidth',2,'displayname','Lattice Struts')
axis equal; hold on
plotV(Vt(fixed,:),'b.','markersize',20,'displayname','Fixed Nodes')
plotV(Vt(forced,:),'r.','markersize',20,'displayname','Loaded Nodes')
theta = linspace(0,2*pi,100);
plot(c_soft(1)+r_soft*sin(theta),c_soft(2)+r_soft*cos(theta),'b--',...
    'linewidth',2,'displayname','Softened')
set(gcf,'position',[360 175 973 774])
set(gca,'FontName','georgia','fontsize',16)
legend('location','northwest')
cmap = flipud(bone(100));
cmap(1:20,:) = [];
colormap(cmap)
title('Nose Cone Geometry, Material, and Loading')
quiver(Xqv,Yqv,Uqv,Vqv,0,'r','linewidth',1.5,'showarrowhead','on',...
    'displayname','Pressure Load')
c = colorbar('location','southoutside');
c.Title.String = 'Normalized Lattice Stiffness';
set(gca,'XColor','none','YColor','none')

%% Now Analyze the Locally Softened Nose Cone
%let's use a linear direct stiffness method

E = Evec_soft';
G = E./(2*(1+nu));

D=struct('m',m,'n',n,'Coord',Coord','Con',Con','Re',Re',...
    'Load',Load','w',w','E',E','G',G','A',A','Iz',Iz','Iy',Iy',...
    'J',J','St',St','be',be');

%Use MSA to solve for displacements
[Q,V,R]=MSA(D);

%extract displacements
U = V(1:2,:)';

%plot undeformed and deformed results
figure
patch('Faces',Ft,'Vertices',Vt,'facecolor','none','edgecolor',[0.8 0.8 0.8],...
    'linewidth',1.5,'displayname','Lattice Struts','handlevisibility','off')
axis equal; hold on
patch('Faces',Ft,'Vertices',Vt+U*10,'facecolor','none','edgecolor','interp',...
    'CData',vecnorm(U*1e3,2,2),'linewidth',2,'displayname','Lattice Struts')
axis equal; hold on
set(gcf,'position',[360 175 973 774])
set(gca,'FontName','georgia','fontsize',16)
set(gca,'XColor','none','YColor','none')
colormap(flipud(brewermap(20,'RdYlBu')))
c = colorbar('location','southoutside');
c.Title.String = 'Displacement [mm]';
title('Uniform Stiffness Cone (disp. magnified 10x)','FontSize',18)
