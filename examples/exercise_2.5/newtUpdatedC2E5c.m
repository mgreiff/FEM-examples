%% Solves the in Fig. 2.6 problem using the total Lagrangian formulation
a = 9;
b = 4;
c = 0;
E = 1;
A = 10;

%%%%%%%%%% Options %%%%%%%%%%%
useImperfectGeometry = 0;
useNonLinearMaterialModel = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmax = 60;
imax = 50;
LIMIT = 1e-4;

edof = [1,1,2,3,4,5,6;
        2,4,5,6,7,8,9];
ex = [0,b;
      b,2*b];
ey = [0,c;
      c,0];
ez = [0,a;
      a,0];

nelm = 2;
ndof = 9;

df = zeros(ndof,1);
df(6) = -0.1;

u = zeros(ndof,1);
f = zeros(ndof,1);

bc = [1,2,3,4,5,7,8,9;
      0,0,0,0,0,0,0,0]';
        
xval = [];
fval = [];

% Re-formulate coordinates
allcord = [0,b,2*b;
           0,c,0;
           0,a,0];

x = reshape(allcord,[9,1]);

xn = x; %Original coordinates
esn = zeros(nelm,1); %Original strain

for n = 1:nmax
    
    f = f + df;
    deltau = zeros(ndof,1);
    strains  = esn;
    
    for i = 1:imax
        
        % Computes the global stiffness matrix
        K = zeros(ndof);
        for ii = 1:nelm
            ind = edof(ii,2:7);
            ec = [ex(ii,:);ey(ii,:);ez(ii,:)];
            l0 = sqrt((ec(:,2) - ec(:,1))'*(ec(:,2) - ec(:,1)));
            Ke = bar3geUpdated( xn(ind), l0, strains(ii), E, A);
            K(ind,ind) = K(ind,ind) + Ke;
        end
        
        % Computes internal forces for bars
        fint = zeros(ndof,1);
        for ii = 1:nelm
            ind = edof(ii,2:7);
            ec = [ex(ii,:);ey(ii,:);ez(ii,:)];
            l0 = sqrt((ec(:,2) - ec(:,1))'*(ec(:,2) - ec(:,1)));
            N = E*A*strains(ii);
            finte = bar3gfUpdated( xn(ind) , l0 , N);
            fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
        end
 
        for ii =1:nelm
            ind = edof(ii,2:7);
            ec = [ex(ii,:);ey(ii,:);ez(ii,:)];
            l0 = sqrt((ec(:,2) - ec(:,1))'*(ec(:,2) - ec(:,1)));
            [ greenStrain, ~ ] = bar3gsUpdated(xn(ind), deltau(ind), l0, E, A);
            strains(ii) = esn(ii) + greenStrain;
        end
        
        r = f - fint;
        du = solveq(K , r , bc);
        deltau = deltau + du;
        x = xn + deltau;

        if norm(r(6)) < LIMIT * norm(df)
            break
        end
        if i == imax
            disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
        end
    end
    
    xn = x;
    esn = strains;
    
    xval = [xval, x];
    fval = [fval, f];
end
%% Plotting
color = 'b';
if useImperfectGeometry
    color = 'r';
end
if useNonLinearMaterialModel
    color = 'g';
end
%% Plots the force-displacement curve
if 1
    figure(1)
    hold on;
    plot(xval(6,:),fval(6,:),[color,'*'])
    xlabel('Elongation, u')
    ylabel('Applied force, f')
    hold off;
end