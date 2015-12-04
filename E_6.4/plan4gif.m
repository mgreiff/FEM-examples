function [ ef ] = plan3gf( ec, t, ed , es )
% Calculating the internal force vector
x = ec(1,:);
y = ec(2,:);

% Compute area of element
eA = 0.5*abs((x(1)*y(2)-y(1)*x(2))+...
             (x(2)*y(3)-y(2)*x(3))+...
             (x(3)*y(4)-y(3)*x(4))+...
             (x(4)*y(1)-y(4)*x(1)));

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

B = B0 + Au * H;
ef = (B'* es) .* (eA * t);
end