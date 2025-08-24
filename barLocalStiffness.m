function k = barLocalStiffness(V)

%V is a 2xdim matrix of the [X Y (Z)] coords of a pair of points

DV = diff(V);           %subtract the endpoints
L = sqrt(sum(DV.^2));   %compute the length of this member
C = DV/L;               %compute the cosine angles Cx Cy Cz
L = C'*C;               %stiffness submatrix lambda
k = [L -L; -L L];       %full stiffness matrix

end