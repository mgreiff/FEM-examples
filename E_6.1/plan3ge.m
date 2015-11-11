function [ Ke ] = plan3ge( ec , t , D , ed , es )
% Calculates the stiffness matrix
x = ec(1,:).*1e-3;
y = ec(2,:).*1e-3;
ed = ed.*1e-3;

% Compute area of element
eA = 0.5*abs(x(1)*(y(2)-y(3))+x(2)*(y(3)-y(1))+x(3)*(y(1)-y(2)));

% 4x6 matrix
H = (1/(2*eA)).*[y(2) - y(3) , 0 , y(3) - y(1) , 0 , y(1) - y(2) , 0;
                 x(3) - x(2) , 0 , x(1) - x(3) , 0 , x(2) - x(1) , 0;
                 0 , y(2) - y(3) , 0 , y(3) - y(1) , 0 , y(1) - y(2);
                 0 , x(3) - x(2) , 0 , x(1) - x(3) , 0 , x(2) - x(1)];

% 3x6 matrix
B0 = [y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0;
      0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1);
      x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),x(2)-x(1)];

% 3x4 matrix
a = H * ed;
Au = [a(1),0,a(3),0;
      0,a(2),0,a(4);
      a(1),a(2),a(3),a(4)];

% 2x2 matrix
S = [es(1) , es(3); es(3) , es(2)];

% 4x4 matrix 
R = [S,zeros(2);zeros(2),S];

K0 = B0'*D*B0;
Ksigma = H'*R*H;
Ku = B0'*D*Au*H+H'*Au'*D*B0+H'*Au'*D*Au*H;
Ke = (K0 + Ku + Ksigma) * (eA * t);
end

