function [ ef ] = plan4gif( ec, t, ed , es )
% Calculating the internal force vector
x = ec(1,:) + ed(1:2:7);
y = ec(2,:) + ed(2:2:8);

% Computes the Jacobian
J11 = @(n) sum(x.*[(n-1),-(n-1),(n+1),-(n+1)]./4);
J12 = @(n) sum(y.*[(n-1),-(n-1),(n+1),-(n+1)]./4);
J21 = @(xi) sum(x.*[(xi-1),-(xi+1),(xi+1),-(xi-1)]./4);
J22 = @(xi) sum(y.*[(xi-1),-(xi+1),(xi+1),-(xi-1)]./4);
J = @(n,xi) [J11(n)  J12(n);
             J21(xi) J22(xi)];

B = @(n, xi) [(n-1),    0,  -(n-1),    0,   (n+1)  ,  0,  -(n+1),   0    ;
                0,   (xi-1),  0,    -(xi+1),   0,  (xi+1),   0,   -(xi-1);
              (xi-1),(n-1) ,-(xi+1),-(n-1), (xi+1),(n+1), -(xi-1),-(n+1)];


F = @(n,xi) B(n,xi)'*es.*det(J(n,xi));

% Gauss points - four point rule OBS ett S per gauss-punkt
p = [-sqrt(3/5),0,sqrt(3/5)];
w = [5/9,8/9,5/9];

ef = 0;
for ii  = 1:length(p)
    for jj = 1:length(p)
        ef = ef + w(ii).*w(jj).*F(p(ii),p(jj)).*t; % dimension 1?
    end
end
end

  