function [] = rhoPlot(DT,rho)

NC = DT.Points;                     %Nodal Coordinates (NC); n_point x 3
LI = edges(DT);                     %List of Edges

% Plot Results
Contours = 24;
C = brewermap(Contours,'Greys');
colormap(C);

%form a matrix of vertices which are unique to each edge
V = [NC(LI(:,1),:); NC(LI(:,2),:)];

%form a matrix of 1D patches by connecting these vertices
F = [(1:length(LI))' (1:length(LI))'+length(LI)];

%form a color data matrix from the element densities
C = [rho;rho];

%plot using patch
patch('faces',F,'vertices',V,'cdata',C,'edgecolor','flat','linewidth',2); hold on
colorbar; axis equal

end
