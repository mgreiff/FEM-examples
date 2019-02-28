function [ S ] = stress1D( E , eG )
%% A non-linear material model
S = (E/2).*log(2.*eG + 1)/sqrt(2.*eG + 1);
end

