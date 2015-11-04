function [ epsilon, N ] = bar3gs( ec , ep , ed )
%Computes (1) Green's stain and (2) the corresponding normal force
x0 = [ec(:,1)',ec(:,2)'];
l0 = sqrt(sum((ec(:,2) - ec(:,1)).^2));
u = ed;
EA = ep(1) .* ep(2);

epsilon = (1/(2.*l0^2)) .* (2*x0 + u) *...
          [eye(3), -eye(3); -eye(3), eye(3)] * u';
N = EA .* epsilon;
end