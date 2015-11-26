function [ S ] = stresscal( eff , E, v)

    F = [eff(1:2), eff(2:3)];

    C = F'*F;
    J = det(F);
    G = E/(2*(1 + v));
    K = E/(3*(1-2*v));

    Cpp = 1; % What is Cpp?
    deltaij = 1; % What is deltaij?

    S = zeros(2);
    for ii = 1:2
        for jj = 1:2
            S(ii,jj) = (K/2).*(J.^2 - 1)/C(ii,jj) +...
                       G.*(J.^(2/3)).*(delta - (Cpp/3)/C(ii,jj));
        end
    end
end

