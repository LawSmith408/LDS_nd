function setupLDS()
%SETUPLDS Add the LDS_nd source, demo, and mesh folders to the MATLAB path.
%
%   Run this once per MATLAB session before executing anything in demos/.
%   It makes the demos independent of the current working directory, so
%   meshes can be loaded by bare filename, e.g.
%
%       setupLDS
%       load small_bunny_mesh2D_vfine.mat
%
%   External dependencies are NOT added here -- see README.md.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'));
addpath(fullfile(root,'demos'));
addpath(fullfile(root,'input_meshes'));

end
