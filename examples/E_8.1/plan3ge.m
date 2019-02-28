function [ Ke ] = plan3ge( ec , t , D , ed , es )
% Calculates the stiffness matrix
x = ec(1,:)';
y = ec(2,:)';

% Compute area of element
eA = 0.5*abs(x(1)*(y(2)-y(3))+x(2)*(y(3)-y(1))+x(3)*(y(1)-y(2)));

% 4x6 matrix
H = (1/(2*eA)).*[y(2) - y(3) , 0 , y(3) - y(1) , 0 , y(1) - y(2) , 0;
                 x(3) - x(2) , 0 , x(1) - x(3) , 0 , x(2) - x(1) , 0;
                 0 , y(2) - y(3) , 0 , y(3) - y(1) , 0 , y(1) - y(2);
                 0 , x(3) - x(2) , 0 , x(1) - x(3) , 0 , x(2) - x(1)];

% 3x6 matrix
B0 = (1/(2*eA)).*[y(2)-y(3),0,          y(3)-y(1),0,            y(1)-y(2),0;
                  0,x(3)-x(2),          0,x(1)-x(3),            0,x(2)-x(1);
                  x(3)-x(2),y(2)-y(3),  x(1)-x(3),y(3)-y(1),    x(2)-x(1),y(1)-y(2)];

% 3x4 matrix
a = H * ed;
Au = [a(1), 0  ,a(3), 0  ;
       0  ,a(2), 0  ,a(4);
      a(2),a(1),a(4),a(3)];

% 2x2 matrix
S = [es(1) , es(3); es(3) , es(2)];

% 4x4 matrix 
R = [S,zeros(2);zeros(2),S];

B = Au * H + B0;

K0 = B'*D*B;
Ksigma = H'*R*H;
Ke = (K0 + Ksigma) .* (eA * t);
end

