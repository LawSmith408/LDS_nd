function U = LDS_Beam_Solver(DT,varargin)
%% A Direct Stiffness method solver for BEAM Lattices
%Lawrence Smith | lasm4254@colorado.edu

%Dependencies:
%https://www.mathworks.com/matlabcentral/fileexchange/27012-matrix-structural-analysis)

%INPUTS:

% DT: a struct containing a triangular or tetrahedral traingulation

%Optional Inputs:
% Evec: [ex1] vector containing the elastic modulus of each bar in the
% mesh where e is the number of edges in the triangulation. If omitted,
% Evec = ones(e,1)

% d_vec: [ex1] vector containing the diameter of each bar in the
% mesh where e is the number of edges in the triangulation. If omitted,
% d_vec = ones(e,1)

%fixed: vector containing the indices of nodes which should be fixed.
%If omitted, all the nodes with the minimum x-coordinate are fixed

%forced: vector containing the indices of the nodes which should be loaded.
%If omited, all nodes with max x coordinate are loaded

%load: a [dx1] vector containing the total force applied to the loaded
%nodes. If omitted, d = [0 -1]' or [0 0 -1]' a unit force in the -vertical direction.
%if dimension is [dxn] where n is the number of loaded nodes, individual
%loads are applied to each n node

NC = DT.Points;         %Nodal Coordinates (NC); n_point x 3
LI = edges(DT);         %List of Edges
dim = size(NC,2);       %dimensionality of probem (2 or 3) 
m = size(LI,1);         %number of edges
n = size(NC,1);         %number of nodes
nu = 0.33;              %Poisson's Ratio

%initialize materials
E = ones(size(LI,1),1);
d_vec = ones(size(LI,1),1);

%initialize load
load = zeros(dim,1);
if dim == 2
    load(end) = -1;
end

% Change these to indicate which points are forced and which are fixed
fixed =  find(NC(:,1)==min(NC(:,1)));       %all nodes where x-coord
forced = find(NC(:,1)==max(NC(:,1)));       %all nodes where x-coord

%parse inputs
if length(varargin)>0
    if ~isempty(varargin{1})
        if numel(varargin{1}) == 1
            E = varargin{1}*ones(1,m);
        else
            E = varargin{1};
        end
    end
end

if length(varargin)>1
    if ~isempty(varargin{2})
        if numel(varargin{2}) == 1
            d_vec = varargin{2}*ones(1,m);
        else
            d_vec = varargin{2};
        end
    end
end

if length(varargin)>2
    if ~isempty(varargin{3})
        fixed = varargin{3};
    end
end

if length(varargin)>3
    if ~isempty(varargin{4})
        forced = varargin{4};
    end
end

if length(varargin)>4
    if ~isempty(varargin{5})
        load = varargin{5};
    end
end

% If this is a 2D problem, add a third column of zeros to the coordinates
if dim==3
    Coord = NC;
elseif dim == 2
    Coord = [NC zeros(size(NC,1),1)]; %coordinates of nodes
end

Con = [LI ones(size(LI,1),2)];    %connectivity of each node and whether rotational dofs are held

%apply fixed boundary conditions
Re = zeros(n,6);                 
Re(fixed,:) = 1;

%apply nodal loads
Load = zeros(n,6); 
Load(forced,1:length(load)) = load(:)'; %Apply load
w = zeros(m,3);

%compute section properties
G = E./(2*(1+nu));
A = pi*d_vec.^2/4;
Iz = pi*d_vec^4/64;
Iy = pi*d_vec^4/64;
J = pi*d_vec^4/32;
St = zeros(n,6);
be = zeros(1,m);

%prepare input structure
D=struct('m',m,'n',n,'Coord',Coord','Con',Con','Re',Re',...
    'Load',Load','w',w','E',E','G',G','A',A','Iz',Iz','Iy',Iy',...
    'J',J','St',St','be',be');

%solve
[~,V,~]=MSA(D);

%extract
if dim==2
U = V(1:2,:)';
elseif dim == 3
U = V(1:3,:)';
end

end