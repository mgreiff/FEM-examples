function [ D ] = dmat1D( E , eG )
%% Computes tangent stiffness
D = E/4*(2*1/(2*eG+1)*1/(sqrt(2*eG+1))+2*log(2*eG+1)*(-1/2)*(2*eG+1)^(-3/2));
end

