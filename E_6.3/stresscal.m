function [ out ] = stresscal( eff , E, v )

    Fmin = [eff(1),eff(2);
            eff(3),eff(4)];
    F = zeros(3);
    F(1:2,1:2) = Fmin;
    F(3,3) = 1; %Material thisckness

    C = F'*F;
    invC = inv(C);
    J = det(F);
    G = E/(2*(1 + v));
    K = E/(3*(1-2*v));

    Cpp = trace(C);

    S = zeros(2);
    for ii = 1:2
        for jj = 1:2
            if ii == jj
                deltaij = 1;
            else
                deltaij = 0;
            end
            S(ii,jj) = (K/2).*(J.^2 - 1).*invC(ii,jj) +...
                       G.*(J.^(-2/3)).*(deltaij - (Cpp/3).*invC(ii,jj));
        end
    end
    out = [S(1,1), S(2,2), S(1,2)]';
end