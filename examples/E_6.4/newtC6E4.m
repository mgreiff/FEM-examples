%%%% Newton rhapson solver for exercise 6.4 %%%%
%% Define geometry and plot initial configuration
load('geom7e4.mat')
if 0
    hold on;
    figure(1)
    for ii = 1:length(ex(:,1))
        ind = [1,2,3,4,1];
        for jj = 1:4
            plot([ex(ii,ind(jj)) , ex(ii, ind(jj + 1))],...
                 [ey(ii,ind(jj)) , ey(ii, ind(jj + 1))], 'r')
        end
    end
end

% Defines additional material properties
E = 210*1e9;
v = 0.3;

%% Sets iteration constants and initial variables

nelm = length(edof(:,1));
ndof = max(max(edof));

nmax = 10;
imax = 20;
LIMIT = 1e-4;

u = zeros(ndof, 1);
fext = zeros(ndof,1);

t = 1;

uval = [];
fval = [];

bc(end-8:end,2) = 2;

% Computes initial fint
fint = zeros(ndof,1);
for ii = 1:nelm
    ec = [ex(ii,:)', ey(ii,:)']';
    ed = u(edof(ii,2:end))';
    [ee , ef] = plan4gis( ec , ed);
    es = cell(4,1);
    for jj = 1:4
        D = mstiff( ef{jj} , E, v );
        es{jj} = D*ee{jj};
    end
    finte = plan4gif( ec, t, ed , es );
    fint(edof(ii,2:end)) = fint(edof(ii,2:end)) + finte;
end
        
for n = 1:nmax

    disp(n)
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
            ed = u(edof(ii,2:end))';
            [ee , ef] = plan4gis( ec , ed);
            D = cell(4,1);
            es = cell(4,1);
            for jj = 1:4
                D{jj} = mstiff( ef{jj} , E, v );
                es{jj} = D{jj}*ee{jj};
            end
            Ke = plan4gie( ec , t , D , ed , es );
            K(edof(ii,2:end),edof(ii,2:end)) = K(edof(ii,2:end),edof(ii,2:end)) + Ke;
        end

        r = fext - fint;
        r(bc(:,1)) = 0;
        du = solveq(K , r , deltaBC);
        u = u + du;

        deltaBC(:,2) = 0;

        fint = zeros(ndof,1);
        for ii = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)']';
            ed = u(edof(ii,2:end))';
            [ee , ef] = plan4gis( ec , ed);
            es = cell(4,1);
            for jj = 1:4
                D = mstiff( ef{jj} , E, v );
                es{jj} = D*ee{jj};
            end
            finte = plan4gif( ec, t, ed , es );
            fint(edof(ii,2:end)) = fint(edof(ii,2:end)) + finte;
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
    ind = [1,2,3,4,1];
    for ii = 1:nelm
        for jj = 1:4
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
    ind = [1,2,3,4,1];
    figure(1)
    hold on;
    for ii = 1:nelm
        defx = ex(ii,:) + u(edof(ii, 2:2:8))';
        defy = ey(ii,:) + u(edof(ii, 3:2:9))';
        for jj = 1:4
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