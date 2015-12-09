%%%% Newton rhapson solver for exercise 6.3 %%%%
%% Define geometry
run('dataC6E3.m')

% Plots the original geometry
if 1
    hold on;
    figure(1)
    for ii = 1:length(ex(:,1))
        ind = [1,2,3,1];
        for jj = 1:3
            plot([ex(ii,ind(jj)) , ex(ii, ind(jj + 1))],...
                 [ey(ii,ind(jj)) , ey(ii, ind(jj + 1))], 'b')
        end
    end
end

%% Sets iteration constants and initial variables
nmax = 30;
imax = 20;
LIMIT = 1e-4;

u = zeros(ndof, 1);
fext = zeros(ndof,1);
df = zeros(ndof,1);
df([3 5]) = 10;

bc = [1 2 6 7;
      0 0 0 0]';

uval = [];
fval = [];

fint = zeros(ndof,1);
for jj = 1:nelm
    ec = [ex(ii,:)', ey(ii,:)']';
    [~, eff] = plan3gs( ec , u(edof(ii,2:7)));
    es = stresscal( eff , E, v );
    finte = plan3gf( ec, t, ed , es' );
    fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
end
        
%% Newton iteration
for n = 1:nmax
    %Applies force
    fext = fext + df;
    disp(fext)
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
            [ ~ , eff ] = plan3gs(ec , ed);
            D = mstiff( eff , E, v );
            es = stresscal( eff , E, v );
            Ke = plan3ge( ec , t , D , ed , es );
            K(edof(ii,2:7),edof(ii,2:7)) = K(edof(ii,2:7),edof(ii,2:7)) + Ke;
        end
        
        r = fext - fint;
        du = solveq(K , r , bc);
        u = u + du;
        
        fint = zeros(ndof,1);
        for jj = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)']';
            [~, eff] = plan3gs( ec , u(edof(ii,2:7)));
            es = stresscal( eff , E, v );
            finte = plan3gf( ec, t, ed , es' );
            fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
        end

        r = fext - fint;
        r(bc(:,1)) = 0;
        
        if i == imax
            disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
        end
    end
    uval = [uval, u];
    fval = [fval, fext];
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
        end
    end
 end