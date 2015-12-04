function [ Ke ] = bar3geUpdated( xn , l0, ee , E, A)
%% Computes the siffness matrix in a 3D bar element
% RETURNS
%     Ke: The stiffness matrix (6x6 array)

% Note that we're using nx and not xtilde
x = xn(4:6) - xn(1:3);

% Linear stiffness matrix
K0 = (E*A/l0.^3) .* [x*x', -x*x'; -x*x', x*x'];

% Initial stress matrix
Ksigma = (E*A*ee/l0)* [eye(3), -eye(3); -eye(3), eye(3)];

Ke = K0 + Ksigma;
end

