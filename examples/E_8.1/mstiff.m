function [ D ] = mstiff( eff , E, v )

    Fmin = [eff(1),eff(2);
            eff(3),eff(4)];
    F = zeros(3);
    F(1:2,1:2) = Fmin;
    F(3,3) = 1; %Material thickness

    C = F'*F;
    invC = inv(C);
    J = det(F);

    G = E/(2*(1 + v));
    K = E/(3*(1-2*v));
    
    Cpp = trace(C);
    
    a1 = K*J.^2 + ((2.*G)/9)*J.^(-2/3)*Cpp;
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
            if i == j
                deltaij = 1;
            else
                deltaij = 0;
            end
            if k == l
                deltakl = 1;
            else
                deltakl = 0;
            end
            D(row, col) = a1.*invC(i,j).*invC(k,l) - ...
                          a2.*(deltaij.*invC(k,l) + deltakl.*invC(i,j)) + ...
                          a3.*(invC(i,k).*invC(j,l) + invC(i,l)*invC(j,k));
        end
    end
end

