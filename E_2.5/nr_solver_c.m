%% Solves the in Fig. 2.6 problem using the total Lagrangian formulation
a = 9;
b = 4;
c = 0;
E = 1;
A = 10; % E*A = 10, and E, A can be set arbitrarily as long as this is acheived

nmax = 30;
imax = 50;
LIMIT = 1e-4;

edof = [1, 1, 2;
        1, 2, 3];
coord = [0,0,0;
         b,c,a;
         2*b,0,0]';

df = zeros(9,1);
df(6) = 0.25;

u = zeros(9,1);
f = zeros(9,1);
x0 = reshape(coord,[9,1]);
x = x0;

nelm = length(edof(:,1));
dof = [1,2,3;4,5,6;7,8,9];

bc = [1,2,3,4,5,7,8,9;
      0,0,0,0,0,0,0,0]';

ep = [E, A];
fval = [];
uval = [];

for n = 1:nmax
    %Computes external force
    f = f - df;
    du = zeros(9,1);
        
    for i = 1:imax
        % Computes the global stiffness matrix
        K = zeros(length(u));
        for ii = 1:nelm
            ec = [coord(:,edof(ii,2)),coord(:,edof(ii,3))];
            indices = (1 + (edof(ii,2)-1)*3:3 + (edof(ii,3)-1)*3);
            ed = u(indices)';
            [~ , es] = bar3gs( ec , ep , ed );
            Ke = bar3ge( ec , ep , ed , es );        
            for jj = 1:2
                for kk = 1:2
                    Kp = Ke((1:1+2)+3*(jj-1), (1:1+2)+3*(kk-1));
                    rows = (1:1+2) + (edof(ii,jj + 1) - 1)*3;
                    cols = (1:1+2) + (edof(ii,kk + 1) - 1)*3;
                    K(rows,cols) = K(rows,cols) + Kp;
                end
            end
        end
        
        fint = zeros(9,1);
        for jj = 1:nelm
            ec = [coord(:,edof(ii,2)),coord(:,edof(ii,3))];
            indices = (1 + (edof(ii,2)-1)*3:3 + (edof(ii,3)-1)*3);
            
            [ee, ~] = bar3gs( ec , ep , u(indices)' );
            
            l0 = sqrt((ec(:,2) - ec(:,1))'*(ec(:,2) - ec(:,1)));
            
            finte = (E * A / l0) * ee * [eye(3), -eye(3); -eye(3), eye(3)] * x(indices);
            fint(indices) = fint(indices) + finte;
        end

        r = f - fint;
        du = solveq(K , r , bc);
        u = u + du;
        x = x0 + u;

        if norm(r(6)) < LIMIT * norm(df)
            break
        end
        if i == imax
            disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
        end
    end
    uval = [uval, u(6)];
    fval = [fval, f(6)];
end
plot(uval,fval,'b*')
xlabel('Elongation, u')
ylabel('Applied force, f')