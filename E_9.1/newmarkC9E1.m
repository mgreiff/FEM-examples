%% Loads and draws undeformed geometry
load('geometryDataC9E1.mat')

% Scales geometry
ex = ex.*1e-3;
ey = ey.*1e-3;
if 1
    hold on;
    figure(1)
    indices = [1,2,3,1];
    for ii = 1:length(edof(:,1))
        for jj = 1:3
            plot([ex(ii,indices(jj)),ex(ii,indices(jj + 1))],...
                 [ey(ii,indices(jj)),ey(ii,indices(jj + 1))], 'b')
        end
    end
    % Plots vertex where pressure is to be applied
    pressureIndex = 64;
    plot(ex(pressureIndex,1),ey(pressureIndex,1),'r*')
    
    % Plots boundary conditions
    plot(ex(1,1),ey(1,1),'kx')
    plot(ex(1,1),ey(1,1),'g+')
    plot(ex(end,2),ey(end,2),'kx')
    plot(ex(end,2),ey(end,2),'g+')
    hold off;
end

%% Sets constants
% Ieration constants (gamma and beta imply trapezoidal rule)
nmax = 80;
imax = 20;
LIMIT = 1e-4;
gamma = 1/2; 
beta = 1/4;
dt = 1;

% Material constants
E = 210e9; % [Pa]
A = 1*1e-6; %    [m^2]
v = 0.3;
t = 1; % [m] Material thickness
rho = 7800; %[kg/m^3]
ep = [E, A];

% Geometry constants
nelm = length(edof(:,1));
ndof = max(max(edof));

%% Newton iteration
% Initial iteration quantities (initial element conditions)
a = zeros(ndof,1);
aval = [];
aN = zeros(ndof,1);
velN = zeros(ndof,1);
accN = zeros(ndof,1);
fpattern = P.*1e6;
fint = zeros(ndof,1);
        
% Begins iteration
for n = 1:nmax
    
    fext = dt.*fpattern; % Applied force increases linearly in time
    
    c1 = 1/(beta.*dt.^2);
    c2 = 1/(beta.*dt);
    c3 = (1-2.*beta)/(2.*beta);
    c4 = gamma/(beta.*dt);
    c5 = (gamma - beta)/beta;
    c6 = dt*(gamma-2*beta)/(2*beta);
    
    d1 = 1;
        
    accDot = c1.*aN + c2.*velN + c3.*accN;
    velDot = c4.*aN + c5.*velN + c6.*accN;
    accPrim = accDot + d1.*velDot;

    Geff = zeros(ndof,1);
    
    for i = 1:imax

        disp([n,i])
        
        %Convergence condition
        if norm(Geff) < LIMIT && i ~= 1
            disp('completed iteration')
            break
        end

        % Computes and assembles the tangent stiffness and mass matrix
        Kt = zeros(ndof);
        M = zeros(ndof);
        fint = zeros(ndof,1);
        for ii = 1:nelm
            ind = edof(ii,2:end);
            ec = [ex(ii,:);ey(ii,:)];
            ed = a(ind);
            [ ~ , eff ] = plan3gs(ec , ed);
            D = mstiff( eff , E, v );
            es = stresscal( eff , E, v );
            Ke = plan3ge( ec , t , D , ed , es );
            Kt(ind, ind) = Kt(ind, ind) + Ke;
            
            Me = plan3gm(ec, t, rho);
            M(ind, ind) = M(ind, ind) + Me;

            finte = plan3gf( ec, t, ed , es );
            fint(ind) = fint(ind) + finte;
        end
        
        % Out of balance forces
        Geff = (c1 + d1*c4).*M*a + fint - fext - M*accPrim;
        
        % Computes displacement increments
        da = solveq(Kt , -Geff , bc);
        Geff(bc(:,1)) = 0;
        
        % Velocity and acceleration
        a = a + da;
        
        acc = c1.*a-accN;
        vel = c4.*a-velN;

    end
    aN = aN + a;
    velN = velN + vel;
    accN = accN + acc;
    
    aval = [aval aN];
end
%% Plot undeformed geometry
 if 1
    ind = [1,2,3,1];
    figure(1)
    hold on;
    for ii = 1:nelm
        defx = ex(ii,:);
        defy = ey(ii,:);
        for jj = 1:3
            plot([defx(ind(jj)) , defx(ind(jj + 1))'],...
                 [defy(ind(jj)) , defy(ind(jj + 1))'], 'b')
        end
    end
 end
%% Plot deformed geometry
 if 1
    ind = [1,2,3,1];
    figure(1)
    hold on;
    index = 80;
    for ii = 1:nelm
        defx = ex(ii,:) + aval(edof(ii, 2:2:6),index)';
        defy = ey(ii,:) + aval(edof(ii, 3:2:7),index)';
        for jj = 1:3
            plot([defx(ind(jj)) , defx(ind(jj + 1))'],...
                 [defy(ind(jj)) , defy(ind(jj + 1))'], 'k')
        end
    end
 end