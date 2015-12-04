% Newton rhapson solver for exercise 6.1
run('data.m')

E = 210*1e9;
v = 0.3;

nelm = length(edof(:,1));
ndof = max(max(edof));

nmax = 5;
imax = 20;
LIMIT = 1e-4;

u = zeros(ndof, 1);
f = zeros(ndof,1);

t = 1;

uval = [];
fval = [];

bc(end-4:end,2) = 10;

if 1
    for n = 1:nmax

        deltaBC = bc;
        for i = 1:imax

            % Creates the global stiffness matrix
            K = zeros(ndof);
            for ii = 1:nelm
                ec = [ex(ii,:)', ey(ii,:)']';
                ed = u(edof(ii,2:7));

                [ ~ , eff ] = plan3gs( ec , ed );
                D = mstiff( eff , E, v );
                es = stresscal( eff , E, v );

                Ke = plan3ge( ec , t , D , ed , es );
                K(edof(ii,2:7),edof(ii,2:7)) = K(edof(ii,2:7),edof(ii,2:7)) + Ke;
            end
        
            fint = zeros(ndof,1);
            for jj = 1:nelm
                ec = [ex(ii,:)', ey(ii,:)']';
                ed = u(edof(ii,2:7));
                
                [ ~ , eff ] = plan3gs( ec , ed );
                es = stresscal( eff , E, v );
                
                finte = plan3gf( ec, t, ed , es );
                fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
            end

            r = -fint;
            du = solveq(K , r , deltaBC);
            u = u + du;
            disp(i)
            disp(u(1:5)')

            deltaBC(:,2) = 0;

            if norm(r) < LIMIT
                break
            end
            if i == imax
                disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
            end
        end
        uval = [uval, u];
        fval = [fval, f];
    end
end

 %% Plots deformed geometry
 if 1
    ind = [1,2,3,1];
    figure(2)
    hold on;
    for ii = 1:nelm
        defx = ex(ii,:) + u(edof(ii, 2:2:6))';
        defy = ey(ii,:) + u(edof(ii, 3:2:7))';
        for jj = 1:3
            plot([defx(ind(jj)) , defx(ind(jj + 1))'],...
                 [defy(ind(jj)) , defy(ind(jj + 1))'], 'b')
            plot(-[defx(ind(jj)) , defx(ind(jj + 1))'],...
                 [defy(ind(jj)) , defy(ind(jj + 1))'], 'b')
            plot([defx(ind(jj)) , defx(ind(jj + 1))'],...
                 -[defy(ind(jj)) , defy(ind(jj + 1))'], 'b')
            plot(-[defx(ind(jj)) , defx(ind(jj + 1))'],...
                 -[defy(ind(jj)) , defy(ind(jj + 1))'], 'b')
        end
    end
 end