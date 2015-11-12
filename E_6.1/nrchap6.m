% Newton rhapson solver for exercise 6.1
clear all; close all;
run('data.m')

% Defines additional material properties
E = 210*1e9; %[Pa]; ?? D-matrisen st?mmer d? med control.mat
v = 0.3;

D = (E/((1+v)*(1-2*v))) .* [1-v, v, 0;   v, 1-v, 0; 0, 0, (1-2*v)/2];
disp(D)

nelm = length(ex(:,1));
ndof = 6 * nelm;

deltaf = 1e-4;
nmax = 20;
imax = 50;
LIMIT = 1e-4;

u = zeros(ndof, 1);
t = 0.1; % OBS! arbitrarily set

for n = 1:nmax
    %Applies force
    %f(6) = f(6) - deltaf;
    
    % Creates the global stiffness matrix
    K = zeros(6 * ndof);
    for ii = 1:nelm
        ec = [ex(ii,:)', ey(ii,:)'];
        ed = u(edof(ii, 2:length(edof(1,:))));
        [~ , es] = plan3ge( ec , t , D , ed , es );
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
    
    % Peforms the equillibrium iteration
    %deltau = zeros(9,1);
    %for i = 1:imax
    %    r = f;
    %    for jj = 1:nelm
    %        ec = coord(:,bars(jj,:));
    %        ed = [deltau(bars(jj,1):bars(jj,1) + 2);...
    %              deltau(bars(jj,2):bars(jj,2) + 2)]';
    %        [ee , ~] = bar3gs( ec , ep , ed );
    %        xtilde = [coord(:,bars(jj,1));coord(:,bars(jj,2))];
    %        l0 = sqrt(sum((ec(:,2) - ec(:,1)).^2));
    %        fadd = (E * A / l0) * ee * [eye(3), -eye(3); -eye(3), eye(3)] * xtilde;
    %        r(bars(jj,1):bars(jj,1) + 2) = r(bars(jj,1):bars(jj,1) + 2) - fadd(1:3);
    %        r(bars(jj,2):bars(jj,2) + 2) = r(bars(jj,2):bars(jj,2) + 2) - fadd(4:6);
    %    end
        
    %    % Solves the system and updates the u-vector
    %    a = solveq(K , f , bc);
    %    deltau = deltau + a;
    %    if norm(r) < LIMIT * norm(f)
    %        break
    %    end
    %    if i == imax
    %        disp(['Warning. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
    %    end
    %end
    %u = u + deltau;
    
    % Plots the applied force as a function in each iteration
    %plot(u(6),f(6),'b*')
end