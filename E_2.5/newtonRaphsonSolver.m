%% Solves the in Fig. 2.6 problem using the total Lagrangian formulation
clear all;
a = 9;
b = 4;
c = 0;
E = 1;
A = 10; % E*A = 10, and E, A can be set arbitrarily as long as this is acheived

% Description of geometry in fig 2.6
bars = [1,2;...
        2,3];
coord = [0,0,0;...
         b,c,a;...
         2*b,0,0];

nelm = length(bars(1,:));

nmax = 20;
imax = 50;
LIMIT = 1e-4;
fval = [];
uval = [];
gval = [];
uold = 0;
f = 0;
for n = 1:nmax
    f = f + deltaf;
    u = uold;
    K = bar3ge(x0, uold, epsilon); % TODO Correct for epsilon in bar3ge
    for i = 1:imax
        epsilon = (1/l0^2) .* (1/2) .* (x0 + xn)' * [eye(), -eye(); eye(), eye()]*u;
        r = 1; % TODO write residual eq
        deltau = K \ r;
        u = u + deltau;
        xn = x0 + u;
        if norm(r) < LIMIT * norm(f)
            break
        end
        if i == imax
            disp('Warning. Reached iteration limit without convergence.')
        end
    end
end