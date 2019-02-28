function [ fint ] = bar3gfUpdated( xn , l0 , N)
%% Computes the internal force vector (see equation 2.19)
% RETURNS
%     f: The internal force vector (6x1 array)
fint = (N/l0) .* [eye(3), -eye(3); -eye(3), eye(3)] * xn;
end