%% Solves the in Fig. 2.6 problem using the total Lagrangian formulation
a = 9;
b = 4;
c = 0.01;
k = 0.3;
E = 1;
A = 10; % E*A = 10, and E, A can be set arbitrarily as long as this is acheived

nmax = 80;
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
ndof = 10;

df = zeros(10,1);
df(6) = -0.1;

u = zeros(10,1);
f = zeros(10,1);

bc = [1,2,3,4,7,8,9,10;
      0,0,0,0,0,0,0,0]';
springInd = [6,10];
        
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
            [~ , es] = bar3gs( ec , ep , ed );
            Ke = bar3ge( ec , ep , ed , es );        
            K(edof(ii,2:7),edof(ii,2:7)) = K(edof(ii,2:7),edof(ii,2:7)) + Ke;
        end
        
        % Adds spring stiffness to K
        Kspring = bar1e( k );
        K(springInd,springInd) = K(springInd,springInd) + Kspring;
        
        % Computes internal forces for bars
        fint = zeros(10,1);
        for jj = 1:nelm
            ec = [ex(ii,:);ey(ii,:);ez(ii,:)];
            [ee, ~] = bar3gs( ec , ep , u(edof(ii,2:7))' );
            
            x = reshape(ec, [6,1]) + u(edof(ii,2:7));
            l0 = sqrt((ec(:,2) - ec(:,1))'*(ec(:,2) - ec(:,1)));
            
            finte = (E * A / l0) * ee * [eye(3), -eye(3); -eye(3), eye(3)] * x;
            fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
        end

        % Computes internal forces for spring
        springFint  = bar1f(k,u(springInd));
        fint(springInd) = fint(springInd) + springFint;
        
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
%% Plots the force-displacement curve
if 1
    figure(1)
    plot(uval(6,:),fval(6,:),'b*')
    xlabel('Elongation, u')
    ylabel('Applied force, f')
end
%% Illustrative movie of the deformation
if 1
    figure(2);
    hold on;
    for ii = 1:nelm
        defCoord = reshape([ex(ii,:);ey(ii,:);ez(ii,:)],[6,1])';
        plot3([defCoord(1),defCoord(4)],...
              [defCoord(2),defCoord(5)],...
              [defCoord(3),defCoord(6)],'b')
    end
    ax = gca;
    ax.NextPlot = 'replaceChildren';

    loops = length(uval(1,:));
    F(loops) = struct('cdata',[],'colormap',[]);
    for ii = 1:loops
        pause(0.1)
        subplot(1,2,1);
        hold on;
        for jj = 1:nelm
            defCoord = reshape([ex(jj,:);ey(jj,:);ez(jj,:)],[6,1]) + uval(edof(jj,2:7),ii);
            plot3([defCoord(1),defCoord(4)],...
                  [defCoord(2),defCoord(5)],...
                  [defCoord(3),defCoord(6)],'b')
        end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis([0,8,-.02,9,-15,9,])
        subplot(1,2,2);
        plot(uval(6,1:ii),fval(6,1:ii),'b*')
        axis([-25,0,-8,0])
        drawnow
        F(ii) = getframe;
    end
end