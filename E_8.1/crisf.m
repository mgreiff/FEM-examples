%clear all; close all;

%% Loads and draws undeformed geometry
load('task1_mesh2.mat')
if 1
    eldrawTri2D(ex, ey, edof)
    title('Undeformed geometry');xlabel('x');ylabel('y');zlabel('z');
end

%% Sets constants
% Ieration constants 
nmax = 20;
imax = 20;
LIMIT = 1e-4;

% Material constants
E = 210e9; % [Pa]
A = 1; %    [mm^2]
v = 0.3;
t = 20; % [m,] Material thickness
ep = [E, A];

% Geometry constants
nelm = length(edof(:,1));
ndof = max(max(edof));

%% Newton iteration
% Initial iteration quantities
a = zeros(ndof,1); 
da = zeros(ndof,1);
dan = zeros(ndof,1);
aold = zeros(ndof,1);
lambda = 0;
psi = 0;
l = 1;

% Loading pattern
P = zeros(ndof,1);
P(64) = 1;

fval = [];
uval = [];

% Sets up initial matrices and constants
S = cell(nelm,0);
D = cell(nelm,0);
for ii = 1:nelm
    ind = edof(ii,2:end);
    ec = [ex(ii,:);ey(ii,:)];
    ed = a(ind);
    [ ~ , eff ] = plan3gs( ec , ed );
    De = mstiff( eff , E, v );
    Se = stresscal( eff , E, v );
    S(ii) = {Se};
    D(ii) = {De};
end

for n = 1:nmax

    G = zeros(ndof,1);

    for i = 1:imax

        %Convergence condition
        if norm(G) < LIMIT && i ~= 1
            disp('completed iteration')
            break
        end

        % Computes and assembles the global tangent stiffness matrix at
        % iteration i.
        % Uses "a", "D" and "S" from the previous iteration, i - 1.
        Kt = zeros(ndof);
        for ii = 1:nelm
            ind = edof(ii,2:end);
            ec = [ex(ii,:);ey(ii,:)];
            ed = a(ind);
            es = S{ii};
            Ke = plan3ge( ec , t , D{ii} , ed , es );
            Kt(ind, ind) = Kt(ind, ind) + Ke;
        end

        % Computes computes pseudo displacement increments
        dag = solveq(Kt , -G , bc);
        dap = solveq(Kt , P , bc);

        da = a - aold;

        % Chooses positive direction when s = 0
        s = sign(dan'*dap);
        if s == 0;
            s = 1;
        end
        
        % Finds a suitable lambda
        increment = 0.01;
        maxCount = 100;
        
        dlambda = s*l/(dap'*dap + psi.*P'*P);
        count = 0;
        while count < maxCount
            count = count + 1;
            
            a1 = dap'*dap + psi.*(P'*P);
            a2 = 2.*dap'*(da + dag) + 2.*psi.*dlambda.*(P'*P);
            a3 = (da + dag)'*(da + dag) + psi.*(dlambda.^2).*(P'*P) - l^2;

            % Solves the equation a1.*lambda^2 + a2*lambda+  a3 = 0
            lambdaSol = roots([a1, a2, a3]);
            if ~isreal(lambdaSol)
                % TODO decrease step length and restart iteration
                dlambda = dlambda.*0.5;
            else
                if lambdaSol(1) ~= lambdaSol(2)
                    % Two solutions to the equation
                    a4 = da'*(da + dag);
                    a5 = da'*dap;
                    cVal = (a4+a5.*lambdaSol)/l.^2;
                    [ ~ , index] = max(cVal);
                    dlambda = lambdaSol(index);
                    disp('lambdaSol(1) ~= lambdaSol(2)')
                    break;
                else
                    dlambda = lambdaSol(1);
                    disp('lambdaSol(1) == lambdaSol(2)')
                    break;
                end
            end
            disp('here')
        end

        % Updates lambda
        lambda = lambda + dlambda;

        % Update displacement vector
        a = a + dag + lambda.*dap;
        
        % Stresses and strains
        S = cell(nelm,0);
        D = cell(nelm,0);
        for ii = 1:nelm
            ind = edof(ii,2:end);
            ec = [ex(ii,:);ey(ii,:)];
            ed = a(ind);
            [ ~ , eff ] = plan3gs( ec , ed );
            De = mstiff( eff , E, v );
            Se = stresscal( eff , E, v );
            S(ii) = {Se};
            D(ii) = {De};
        end
        
        % Compute internal forces
        fint = zeros(ndof,1);
        for ii = 1:nelm
            ind = edof(ii,2:end);
            ec = [ex(ii,:);ey(ii,:)];
            ed = a(ind);
            es = S{ii};
            finte = plan3gf( ec, t, ed , es' );
            fint(ind) = fint(ind) + finte;
        end
 
        G = fint - lambda.*P;
        f = - lambda.*P;
    end
    
    dan = a - aold;
    aold = a;

    uval = [uval, a];
    fval = [fval, f];
end