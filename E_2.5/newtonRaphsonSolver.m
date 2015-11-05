%% Solves the in Fig. 2.6 problem using the total Lagrangian formulation
clear all;
a = 9;
b = 4;
c = 0;
E = 1;
A = 10; % E*A = 10, and E, A can be set arbitrarily as long as this is acheived

deltaf = 1e-4;
nmax = 20;
imax = 50;
LIMIT = 1e-4;

% Description of geometry in fig 2.6
bars = [1,2; 2,3];
coord = [0,0,0; b,c,a; 2*b,0,0]';
deltau = zeros(size(coord));
nelm = length(bars(1,:));
ndof = 3;

fval = [];
uval = [];
uold = zeros(9,1);
f = zeros(9,1);

ep = [E, A];

% Boundary conditions
bc = [1,2,3,4,5,7,8,9;0,0,0,0,0,0,0,0]';
disp(bc)

for n = 1:nmax
    f(6) = f(6) - deltaf; %Applied force as in fig 2.6
    u = uold;
    K = zeros(3 * ndof);

    for ii = 1:nelm
        ec = coord(:,bars(ii,:));
        ed = [deltau(:,bars(ii,1)); deltau(:,bars(ii,2))]';
        [~ , es] = bar3gs( ec , ep , ed );
        Ke = bar3ge( ec , ep , ed , es );
        for jj = 1:2
            for kk = 1:2
                Kp = Ke(ndof*jj - (ndof - 1):ndof*jj, ndof*kk - (ndof - 1):ndof*kk);
                rows = ndof * bars(ii,jj) - (ndof - 1):ndof * bars(ii,jj);
                cols = ndof * bars(ii,kk) - (ndof - 1):ndof * bars(ii,kk);
                K(rows,cols) = K(rows,cols) + Kp;
            end
        end
    end

    for i = 1:imax
        r = f;
        for jj = 1:nelm
            ec = coord(:,bars(jj,:));
            ed = [deltau(:,bars(jj,1)); deltau(:,bars(jj,2))]';
            [ee , ~] = bar3gs( ec , ep , ed );
            xtilde = [coord(:,bars(jj,1));coord(:,bars(jj,2))];
            l0 = sqrt(sum((ec(:,2) - ec(:,1)).^2));
            fadd = (E * A / l0) * ee * [eye(3), -eye(3); -eye(3), eye(3)] * xtilde;
            r(bars(jj,1):bars(jj,1) + 2) = r(bars(jj,1):bars(jj,1) + 2) - fadd(1:3);
            r(bars(jj,2):bars(jj,2) + 2) = r(bars(jj,2):bars(jj,2) + 2) - fadd(4:6);
        end
        a = solveq(K , f , bc);
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