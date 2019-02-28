function [ es, ee ] = bar3gsNonLin( ec , ep , ed )
%% Computes the Green's stain and normal force in a 3D bar element
% RETURNS
%     Eg: The Green's stain in the 3D bar element (float)
%     N: The normal force in the 3D bar element (float)
x0 = ec(:,2) - ec(:,1);
l0 = sqrt(x0'*x0);
u = (ed(4:6) - ed(1:3))';
l = sqrt((x0 + u)'*(x0 + u));
E = ep(1);
A = ep(2);

ee = (l^2-l0^2)/(2*l0^2);

S = stress1D(E, ee);
es = A*S;
end