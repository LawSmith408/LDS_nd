function [S] = deformedStressPlot(DT,D,eMat,plotme)

NC = DT.Points;                     %Nodal Coordinates (NC); n_point x 3
LI = edges(DT);                     %List of Edges

if length(eMat)==1
    eMat = eMat*ones(length(LI),1);
end

%compute element lengths
L0 = arrayfun(@(i) sqrt(sum(diff(NC(LI(i,:),:)).^2)), 1:size(LI,1));
L1 = arrayfun(@(i) sqrt(sum(diff(NC(LI(i,:),:)+D(LI(i,:),:)).^2)), 1:size(LI,1));
dL = L1-L0;                        %change in length
S = (dL(:)./L0(:)).*eMat(:);      %axial stress = strain * modulus

if plotme
% Plot Results
Contours = 24;
Scolor = floor(rescale([-max(abs(S)); S(:); max(abs(S));],1,Contours));
C = brewermap(Contours,'RdBu');
colormap(C);

x = NC(:,1)+D(:,1);
y = NC(:,2)+D(:,2);

for i = 1:length(LI)
   %plot deformed elements
   plot([x(LI(i,1)) x(LI(i,2))],[y(LI(i,1)) y(LI(i,2))],...
       'linewidth',2,'color',C(Scolor(i),:)); hold on
end
axis equal

end

end