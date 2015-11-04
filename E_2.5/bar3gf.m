function [ force ] = bar3gf( ec , ed , es)
%Computes the internal force vector (see equation 2.19)
N = es(1);
l0 = sqrt(sum((ec(:,2) - ec(:,1)).^2));
x0 = [ec(:,1)',ec(:,2)'];
x = x0 + ed;
force = ((N/l0) .* [eye(3), -eye(3); -eye(3), eye(3)] * x')';
end