function [ D ] = mstiff( eff , E, v )

    F = [eff(1:2), eff(2:3)];

    C = F'*F;
    J = det(F);

    G = E/(2*(1 + v));
    K = E/(3*(1-2*v));
    
    Cpp = 1; % What is Cpp?
    deltaij = 1; % What is deltaij?
    deltakl = 1; % What is deltakl?
    
    a1 = K*J.^2 + ((2.*G)/9)*J^(-2/3)*Cpp;
    a2 = ((2.*G)/3)*J^(-2/3);
    a3 = (G/3)*J^(-2/3)*Cpp - (K/2)*(J^2 - 1);
    
    D = zeros(3);
    
    indexPairs = [1,1;2,2;1,2];
    for row = 1:3
        for col = 1:3
            i = indexPairs(row, 1);
            j = indexPairs(row, 2);
            k = indexPairs(col, 1);
            l = indexPairs(col, 2);
            D(row, col) = a1/(C(i,j)*C(k,l)) - ...
                          a2*(deltaij/C(k,l) + deltakl/C(i,j)) + ...
                          a3*(1/(C(i,j)*C(k,l)) + 1/(C(k,l)*C(i,j)));
        end    
    end
end

