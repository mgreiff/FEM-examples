function [ Eg, N ] = bar3gs( ec , ep , ed )
%% Computes the Green's stain and normal force in a 3D bar element
% RETURNS
%     Eg: The Green's stain in the 3D bar element (float)
%     N: The normal force in the 3D bar element (float)
x0 = [ec(:,1)',ec(:,2)'];
l0 = sqrt(sum((ec(:,2) - ec(:,1)).^2));
EA = ep(1) .* ep(2);

Eg = (1/(2.*l0^2)) .* (2*x0 + ed) * [eye(3), -eye(3); -eye(3), eye(3)] * ed';
N = EA .* Eg;
end