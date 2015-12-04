function [ Ke ] = bar3ge( ec , ep , ed , es )
%% Computes the siffness matrix in a 3D bar element
% RETURNS
%     Ke: The stiffness matrix (6x6 array)
x0 = ec(:,2) - ec(:,1);
l0 = sqrt(x0'*x0);
EA = ep(1) * ep(2);
N = es(1);
u = (ed(4:6) - ed(1:3))';

% Linear stiffness matrix
K0 = (EA/l0^3) .* [x0*x0', -x0*x0'; -x0*x0', x0*x0'];

% Initial displacement matrix
Ku = (EA/l0^3) .* [  x0*u' + u*x0' + u*u', -(x0*u' + u*x0' + u*u');...
                   -(x0*u' + u*x0' + u*u'),  x0*u' + u*x0' + u*u'];

% Initial stress matrix
Ksigma = (N/l0)* [eye(3), -eye(3); -eye(3), eye(3)];

Ke = K0 + Ku + Ksigma;
end

