% Newton rhapson solver for exercise 6.1
load('geometryC6E2.mat')

% Defines additional material properties
E = 210*1e9; %[Pa]; 
v = 0.3;

nelm = length(edof(:,1));
ndof = max(max(edof));

nmax = 10;
imax = 40;
LIMIT = 1e-2;

u = zeros(ndof, 1);
fext = zeros(ndof,1);

t = 1;

uval = [];
fval = [];

bc(end-4:end,2) = 2;

% Computes initial fint
fint = zeros(ndof,1);
for ii = 1:nelm
    ec = [ex(ii,:)', ey(ii,:)']';
    ed = u(edof(ii,2:7));
    [ee , eff] = plan3gs( ec , u(edof(ii,2:7)));
    D = mstiff( eff , E, v );
    es = D*ee;
    finte = plan3gf( ec, t, ed , es );
    fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
end
        
for n = 1:nmax

    deltaBC = bc;
    r = 1;
    
    for i = 1:imax
        
        if norm(r) < LIMIT
            break
        end
        
        % Creates the global stiffness matrix
        K = zeros(ndof);
        for ii = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)']';
            ed = u(edof(ii,2:7));
            [ ee , eff ] = plan3gs(ec , ed);
            D = mstiff( eff , E, v );
            es = stresscal( eff , E, v );
            Ke = plan3ge( ec , t , D , ed , es );
            K(edof(ii,2:7),edof(ii,2:7)) = K(edof(ii,2:7),edof(ii,2:7)) + Ke;
        end

        r = fext - fint;
        r(bc(:,1)) = 0;
        du = solveq(K , r , deltaBC);
        u = u + du;

        deltaBC(:,2) = 0;

        fint = zeros(ndof,1);
        for ii = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)']';
            ed = u(edof(ii,2:7));
            %[ee, ~] = plan3gs( ec , u(edof(ii,2:7)));
            es = stresscal( eff , E, v );
            finte = plan3gf( ec, t, ed , es );
            fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
        end
        
        r = fext - fint;
        r(bc(:,1)) = 0;
        
        if i == imax
            disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
        end
    end
    uval = [uval, u];
end
%% Plots the elements
if 1
    figure(1)
    hold on;
    ind = [1,2,3,1];
    for ii = 1:nelm
        for jj = 1:3
            plot([ex(ii,ind(jj)) , ex(ii,ind(jj + 1))'],...
                 [ey(ii,ind(jj)) , ey(ii,ind(jj + 1))'], 'b')
            plot(-[ex(ii,ind(jj)) , ex(ii,ind(jj + 1))'],...
                 [ey(ii,ind(jj)) , ey(ii,ind(jj + 1))'], 'b')
            plot([ex(ii,ind(jj)) , ex(ii,ind(jj + 1))'],...
                 -[ey(ii,ind(jj)) , ey(ii,ind(jj + 1))'], 'b')
            plot(-[ex(ii,ind(jj)) , ex(ii,ind(jj + 1))'],...
                 -[ey(ii,ind(jj)) , ey(ii,ind(jj + 1))'], 'b')
        end
    end
    hold off;
end
 %% Plots deformed geometry
 if 1
    ind = [1,2,3,1];
    figure(1)
    hold on;
    for ii = 1:nelm
        defx = ex(ii,:) + u(edof(ii, 2:2:6))';
        defy = ey(ii,:) + u(edof(ii, 3:2:7))';
        for jj = 1:3
            plot([defx(ind(jj)) , defx(ind(jj + 1))'],...
                 [defy(ind(jj)) , defy(ind(jj + 1))'], 'r')
            plot(-[defx(ind(jj)) , defx(ind(jj + 1))'],...
                 [defy(ind(jj)) , defy(ind(jj + 1))'], 'r')
            plot([defx(ind(jj)) , defx(ind(jj + 1))'],...
                 -[defy(ind(jj)) , defy(ind(jj + 1))'], 'r')
            plot(-[defx(ind(jj)) , defx(ind(jj + 1))'],...
                 -[defy(ind(jj)) , defy(ind(jj + 1))'], 'r')
        end
    end
 end