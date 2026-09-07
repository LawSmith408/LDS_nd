function [k,L] = barLocalStiffness(V)
%BARLOCALSTIFFNESS Local stiffness matrix of a 2D or 3D bar (truss) element.
%
%   k = barLocalStiffness(V) returns the [2*dim x 2*dim] unit-stiffness
%   matrix for the bar whose endpoints are the rows of V, a [2 x dim]
%   matrix of [X Y (Z)] coordinates. The result carries no material or
%   section information: scale it by EA/L to obtain the physical element
%   stiffness.
%
%   [k,L] = barLocalStiffness(V) also returns the element length L.
%
%   DOF ordering within k is node-major: [u1x u1y (u1z) u2x u2y (u2z)].

DV = diff(V);           %vector from endpoint 1 to endpoint 2
L  = sqrt(sum(DV.^2));  %length of this member
C  = DV/L;              %direction cosines Cx Cy (Cz)
Lam = C'*C;             %stiffness submatrix lambda
k  = [Lam -Lam; -Lam Lam];

end
