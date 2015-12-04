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
f = zeros(ndof,1);
df = zeros(ndof,1);
df([3 5]) = -1e8;

bc = [1 2 6 7;
      0 0 0 0]';

uval = [];
fval = [];

%% Newton iteration
for n = 1:nmax
    %Applies force
    f = f + df;
    
    for i = 1:imax
        % Creates the global stiffness matrix
        K = zeros(ndof);
        for ii = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)'];
            ed = u(edof(ii,2:7));
            [ ee , ~ ] = plan3gs(ec' , ed);
            D = (E/((1+v)*(1-2*v))) .* [1-v, v, 0;   v, 1-v, 0; 0, 0, (1-2*v)/2];
            es = D * ee; % Compute second order poila-kirchoff stress ??
            Ke = plan3ge( ec' , t , D , ed , es );
            K(edof(ii,2:7),edof(ii,2:7)) = K(edof(ii,2:7),edof(ii,2:7)) + Ke;
        end
        
        fint = zeros(ndof,1);
        for jj = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)']';
            
            [ee, ~] = plan3gs( ec , u(edof(ii,2:7)));
            A = 0.5*abs(ec(1,1)*(ec(2,2)-ec(2,3))+...
                ec(1,2)*(ec(2,3)-ec(2,1))+...
                ec(1,3)*(ec(2,1)-ec(2,2)));
            es = A.*ee;
            finte = plan3gf( ec, t, ed , es );
            fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
        end

        r = f - fint;
        du = solveq(K , r , bc);
        u = u + du;

        if norm(r([3,5])) < LIMIT * norm(df)
            break
        end
        if i == imax
            disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
        end
    end
    uval = [uval, u];
    fval = [fval, f];
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