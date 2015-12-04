function [ Eg, N ] = bar3gsUpdated( xn , deltau, l0, E, A )
%% Computes the Green's stain and normal force in a 3D bar element
% RETURNS
%     Eg: The Green's stain in the 3D bar element (float)
%     N: The normal force in the 3D bar element (float)
Eg = (1/(2.*l0.^2)) .* (2.*xn + deltau)' * [eye(3), -eye(3); -eye(3), eye(3)] * deltau;
N = E.*A.* Eg;
end