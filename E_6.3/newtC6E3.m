%%%% Newton rhapson solver for exercise 6.3 %%%%
%% Define geometry
ex = [0, 10, 10;
      0, 10,  0];
  
ey = [0, 0,  10;
      0, 10, 10];

edof = [1 1 2 3 4 5 6;
        2 1 2 5 6 7 8];

E = 210*1e9;
v = 0.3;
t = 1;

nelm = length(edof(:,1));
ndof = max(max(edof));

%% Options
CompNonLin = 1;
if CompNonLin
    color = 'g';
else
    color = 'r';
end

%% Sets iteration constants and initial variables
nmax = 100;
imax = 30;
LIMIT = 1e-4;

u = zeros(ndof, 1);
fext = zeros(ndof,1);
df = zeros(ndof,1);
df([3 5]) = -2e9;

bc = [1 2 4 7;
      0 0 0 0]';

uval = [];
fval = [];

%% Newton iteration
for n = 1:nmax
    %Applies force
    fext = fext + df;
    
    for i = 1:imax
        disp([n,i])
        
        % Creates the global stiffness matrix
        K = zeros(ndof);
        fint = zeros(ndof,1);
        for ii = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)']';
            ed = u(edof(ii,2:7));
            [ ee , eff ] = plan3gs(ec , ed);
            if CompNonLin
                D = mstiff( eff , E, v );
                es = stresscal( eff , E, v );
            else
                D = (E/((1+v)*(1-2*v))) .* [1-v, v, 0;   v, 1-v, 0; 0, 0, (1-2*v)/2];
                es = D * ee;
            end
            Ke = plan3ge( ec , t , D , ed , es );
            K(edof(ii,2:7),edof(ii,2:7)) = K(edof(ii,2:7),edof(ii,2:7)) + Ke;

            finte = plan3gf( ec , t , ed , es );
            fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
        end

        r = fext - fint;
        du = solveq(K , r , bc);
        u = u + du;
        
        r(bc(:,1)) = 0;
        if norm(r) < LIMIT * norm(df)
            break
        end
        if i == imax
            disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
        end
    end
    uval = [uval, u];
    fval = [fval, fext];
end
%% Plots the original geometry
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

%% Plots deformed geometries
 if 1
    ind = [1,2,3,1];
    figure(1)
    hold on;
    for ii = 1:nelm
        defx = ex(ii,:) + u(edof(ii, 2:2:6))';
        defy = ey(ii,:) + u(edof(ii, 3:2:7))';
        for jj = 1:3
            plot([defx(ind(jj)) , defx(ind(jj + 1))'],...
                 [defy(ind(jj)) , defy(ind(jj + 1))'], color)
        end
    end
 end
 
%% Plots force displacement
 if 1
    figure(2)
    indices = [3];
    for ii = 1:length(indices)
        plot(uval(indices(ii),:), fval(indices(ii),:),color)
    end
 end