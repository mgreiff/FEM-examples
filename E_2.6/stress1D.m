function [ S ] = stress1D( E , eG )
%% A non-linear material model
S = (E/2).*ln(2.*eG + 1)/sqrt(2.*eG + 1);
end

