%% Solves the in Fig. 2.6 problem using the total Lagrangian formulation
a = 9;
b = 4;
c = 0.01;
k1 = 0.3;
k2 = 0.3;
E = 1;
A = 10;

%%%%%%%%%% Options %%%%%%%%%%%
useImperfectGeometry = 1;
useNonLinearMaterialModel = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmax = 80;
imax =30;
LIMIT = 1e-2;

edof = [1,1,2,3,4,5,6;
        2,4,5,6,7,8,9];
ex = [0,b;
      b,2*b];
ey = [0,c;
      c,0];
ez = [0,a;
      a,0];

% Use imperfect geometry as described in f)
if useImperfectGeometry
    displacement = 1;
    ex = [0,b+displacement;
          b+displacement,2*b];
end

nelm = 2;
ndof = 11;

df = zeros(ndof,1);
df(6) = -0.1;

u = zeros(ndof,1);
f = zeros(ndof,1);

bc = [1,2,3,4,7,8,9,10,11;
      0,0,0,0,0,0,0,0 ,0]';
spring1Ind = [5,10];
spring2Ind = [6,11];

ep = [E, A];
fval = [];
uval = [];

for n = 1:nmax
    
    f = f + df;
    du = zeros(ndof,1);
        
    for i = 1:imax
        % Computes the global stiffness matrix
        K = zeros(ndof);
        for ii = 1:nelm
            ec = [ex(ii,:);ey(ii,:);ez(ii,:)];
            ed = u(edof(ii,2:7))';
            if useNonLinearMaterialModel
                [~ , ee] = bar3gsNonLin( ec , ep , ed );
                Ke = bar3geNonLin( ec , ep , ed , ee );
            else
                [es , ~] = bar3gs( ec , ep , ed );
                Ke = bar3ge( ec , ep , ed , es ); 
            end       
            K(edof(ii,2:7),edof(ii,2:7)) = K(edof(ii,2:7),edof(ii,2:7)) + Ke;
        end
        
        % Adds spring stiffness to K
        K(spring1Ind,spring1Ind) = K(spring1Ind,spring1Ind) + bar1e( k1 );
        K(spring2Ind,spring2Ind) = K(spring2Ind,spring2Ind) + bar1e( k2 );
         
        % Computes internal forces for bars
        fint = zeros(ndof,1);
        for jj = 1:nelm
            ec = [ex(ii,:);ey(ii,:);ez(ii,:)];
            ed = u(edof(ii,2:7))';
            if useNonLinearMaterialModel
                [es, ~] = bar3gsNonLin(ec , ep , ed );
                finte = bar3gf(ec , ed , es);
            else
                [es, ~] = bar3gs(ec, ep, ed);
                finte = bar3gf(ec , ed , es);
            end
            fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
        end

        % Computes internal forces for spring
        fint(spring1Ind) = fint(spring1Ind) + bar1f(k1,u(spring1Ind));
        fint(spring2Ind) = fint(spring2Ind) + bar1f(k2,u(spring2Ind));
        
        r = f - fint;
        du = solveq(K , r , bc);
        u = u + du;

        if norm(r(6)) < LIMIT * norm(df)
            break
        end
        if i == imax
            disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
        end
    end
    uval = [uval, u];
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
    plot(uval(6,:),fval(6,:),[color,'*'])
    xlabel('Elongation, u')
    ylabel('Applied force, f')
    hold off;
end
%% Illustrative movie of the deformation
if 1
    figure(2);
    hold on;
    for ii = 1:nelm
        defCoord = reshape([ex(ii,:);ey(ii,:);ez(ii,:)],[6,1])';
        plot3([defCoord(1),defCoord(4)],...
              [defCoord(2),defCoord(5)],...
              [defCoord(3),defCoord(6)],color)
    end
    ax = gca;
    ax.NextPlot = 'replaceChildren';

    loops = length(uval(1,:));
    F(loops) = struct('cdata',[],'colormap',[]);
    for ii = 1:loops
        pause(0.025)
        subplot(1,2,1);
        hold on;
        for jj = 1:nelm
            defCoord = reshape([ex(jj,:);ey(jj,:);ez(jj,:)],[6,1]) + uval(edof(jj,2:7),ii);
            plot3([defCoord(1),defCoord(4)],...
                  [defCoord(2),defCoord(5)],...
                  [defCoord(3),defCoord(6)],color)
        end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis([0,8,-.02,0.05,-15,9])
        subplot(1,2,2);
        plot(uval(6,1:ii),fval(6,1:ii),color)
        axis([-25,0,df(6)*nmax,0])
        drawnow
        F(ii) = getframe;
    end
end