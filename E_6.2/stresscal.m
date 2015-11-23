function [ S ] = stresscal( eff , ed )
F = [eff(1:2), eff(2:3)];
C = F'*F;
S = eye(2) % To be derived
end

