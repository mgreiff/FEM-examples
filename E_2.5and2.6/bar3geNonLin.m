function [ Ke ] = bar3geNonLin( ec , ep , ed , ee )
%% Computes the siffness matrix in a 3D bar element
% RETURNS
%     Ke: The stiffness matrix (6x6 array)

x0 = ec(:,2) - ec(:,1);
l0 = sqrt(x0'*x0);
E = ep(1);
A = ep(2);
u = (ed(4:6) - ed(1:3))';

x = x0 + u;

D = dmat1D(E, ee);
S = stress1D(E, ee);

Ke = (A/l0).*(D*[x;-x]*[x',-x']/l0.^2 + S*[eye(3),eye(3);eye(3),eye(3)]);
end

