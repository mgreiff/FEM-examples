function [ fint ] = bar3gf( ec , ed , es)
%% Computes the internal force vector (see equation 2.19)
% RETURNS
%     f: The internal force vector (6x1 array)
N = es;
l0 = sqrt(sum((ec(:,2) - ec(:,1)).^2));
x0 = [ec(:,1)',ec(:,2)'];
fint = (N/l0) .* [eye(3), -eye(3); -eye(3), eye(3)] * (x0 + ed)';
end