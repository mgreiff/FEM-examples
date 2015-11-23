clear all; close all;
load('task1_mesh2.mat')

%% Draw undeformed geometry
if 1
    eldrawTri2D(ex, ey, edof)
    title('Undeformed geometry');xlabel('x');ylabel('y');zlabel('z');
end

% Sets ieration constants 
nmax = 20;
imax = 20;
LIMIT = 1e-4;

% Material constants
E = 210e9; % [Pa]
A = 1; %    [mm^2]
v=0.3;

ep = [E, A];

nelm = length(edof(:,1));
ndof = numel(Dof);

bc = [1:ndof;
      zeros(1,21)]';
bc([3,6,9,12],:) = [];

% Declares initial states

df = zeros(ndof,1);
df(3) = -1000;

fval = [];
uval = [];

% Initial iteration quantities
a = zeros(ndof,1); 
lambda = 1; % scalar, set arbitarily?
    
for n = 1:nmax
    
    for i = 1:imax

        % Computes the global stiffness matrix at iteration i
        % Use a, compute and use S

        % Computes computes pseudo displacement increments
        
        % Load parameter
        a1 = dap'*dap + psi.*(P'*P);
        a2 = 2.*dap'*(da + dag) + 2.*psi.*dlambda.*(P'*P);
        a3 = (da + dag)'*(da + dag) + psi.*(dlambda.^2).*(P'*P) - l^2;
        
        % Solves the equation a1.*dlambda^2 + a2*dlambda+  a3 = 0
        lambdaSol = roots([a1, a2, a3]);
        if ~isreal(lambdaSol)
            % TODO decrease step length and restart iteration
        end
        if lambdaSol(1) ~= lambdaSol(2)
            % Two solutions to the equation
            a4 = da'*(da + dag);
            a5 = da'*dap;
            cVal = (a4+a5.*lambdaSol)/l.^2;
            [ ~ , index] = max(cVal);
            dlambda = dlambda +  lambdaSol(index);
        else
            dlambda = dlambda +  lambdaSol(1);
        end
        
        % Update displacement vector
        a = a + dag + dlambda.*dap;
        
        % Stresses and strains
        fint = 

        % Compute internal forces
        G = fint - 
    end
    u = u + du;
    f = f + xi.*df;
    
    uval = [uval, u];
    fval = [fval, f];
end