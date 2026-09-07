%% Vectorized 2D/3D Linear Direct Stiffness Method Simulator on Arbitrary Tets/Tris
% Lawrence Smith | lasm4254@colorado.edu

%INPUTS:

% DT: a struct containing a triangular or tetrahedral traingulation

%Optional Inputs:
% Evec: [ex1] vector containing the elastic modulus*area of each bar in the
% mesh where e is the number of edges in the triangulation. If omitted,
% Evec = ones(e,1)

%fixed: vector containing the indices of nodes which should be fixed.
%If omitted, all the nodes with the minimum x-coordinate are fixed

%forced: vector containing the indices of the nodes which should be loaded.
%If omited, all nodes with max x coordinate are loaded

%loadVec: a [dx1] vector containing the total force applied to the loaded
%nodes. If omitted, d = [0 -1]' or [0 0 -1]' a unit force in the -vertical direction.
%if dimension is [dxn] where n is the number of loaded nodes, individual
%loads are applied to each n node

%OUTPUTS:

%D: A matrix containing the displacements at all the nodes
%C: A scalar compliance measure of the structure


function [D,C] = LDS_Bar_Solver(DT,varargin)

NC = DT.Points;                     %Nodal Coordinates (NC); n_point x 3
LI = edges(DT);                     %List of Edges
dim = size(NC,2);                   %dimensionality of probem (2 or 3) 
nel_k = (2*dim)^2;                  %No. of Elements in local stiffness matri
ID = reshape(1:numel(NC),dim,[])';  %Indexing schema

%initialize materials
eMat = ones(size(LI,1),1);

%initialize load vector
loadVec = zeros(dim,1);
if dim == 2
    loadVec(end) = -5;
else
    loadVec(end) = -100;
end

% Change these to indicate which points are forced and which are fixed
fixed =  find(NC(:,1)==min(NC(:,1)));       %all nodes where x-coord
forced = find(NC(:,1)==max(NC(:,1)));       %all nodes where x-coord

% Parse optional inputs
if length(varargin)>0
    if ~isempty(varargin{1})
        if numel(varargin{1}) == 1
            eMat = 0*eMat  + 1*varargin{1};
        else
            eMat = varargin{1};
        end
    end
end

if length(varargin)>1
    if ~isempty(varargin{2})
        fixed = varargin{2};
    end
end

if length(varargin)>2
    if ~isempty(varargin{3})
        forced = varargin{3};
    end
end

if length(varargin)>3
    if ~isempty(varargin{4})
        loadVec = varargin{4};
    end
end

% By process of elimination determine the non-fixed DOF
free  =  ID';
free(:,fixed) = [];
free = free(:);

% Construct load vector
F = zeros(size(NC))';
if size(loadVec,2)==1 %if load direction is constant across all loaded nodes
    F(:,forced) = F(:,forced)+loadVec(:)./numel(forced);  %distribute the total load over the loaded nodes
else
    F(:,forced) = loadVec; %per-node loads supplied directly
end
F = F(:);        %reshape into a column vector

IND = zeros(size(LI,1)*nel_k,2);    % index matr used for local stiffness
K1 = zeros(size(LI,1)*nel_k,1);     %empty accumulation matrix

% compute all element stiffness matrices
tic
for j = 1:size(LI,1)
    [j1, j2] = meshgrid(reshape(ID(LI(j,:),:)',[],1)');
    IND((j-1)*nel_k+1:j*nel_k,:) = [j1(:) j2(:)];
    V = NC(LI(j,:),:);
    k = eMat(j)*barLocalStiffness(V);
    K1((j-1)*nel_k+1:j*nel_k) = k(:);
end

% Assemble global stiffness matrix as a sparse matrix
K = sparse(IND(:,1),IND(:,2),K1,numel(NC),numel(NC));

% Solve the reduced system. Preconditioned conjugate gradients is usually
% faster than backslash, but it can stagnate on badly conditioned systems
% (e.g. a SIMP design with a wide Emin:Emax ratio), so fall back to a
% direct solve whenever it fails to converge.
Kff = K(free,free);
[d,flag,relres,iter] = pcg(Kff,F(free),1e-8,1e4);
if flag ~= 0
    warning('LDS:pcgFailed',...
        ['pcg did not converge (flag %d, relres %.2e after %d iters); '...
         'falling back to a direct solve.'],flag,relres,iter);
    d = Kff\F(free);
end

% Augment displacement matrix with known DOF
D = zeros(numel(NC),1);
D(free) = d;                %Assemble full Displacement Matrix
C = F'*D;                   %compute compliance
D = reshape(D,dim,[])';     %reshape into same dimensionality as NC

fprintf('Solved %.i DOF in %.3fs, Compliance = %.2e\n',numel(d),toc,C)

end
