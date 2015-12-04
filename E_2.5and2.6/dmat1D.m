function [ D ] = dmat1D( E , eG )
%% Computes tangent stiffness
D = (E/2).*log(2.*eG + 1)/sqrt(2.*eG + 1);
end

