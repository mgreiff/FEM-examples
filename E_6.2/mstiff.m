function [ D ] = mstiff( eff , ed )
F = [eff(1:2), eff(2:3)];
C = F'*F;
D = eye(2) % To be derived
end

