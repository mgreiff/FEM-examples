function [ ee , ef ] = plan3gs( ec , ed )
% Calculates Green?s strains and the deformation gradient
% Calculating the internal force vector
x = ec(1,:);
y = ec(2,:);

% Compute area of element
A = 0.5 * abs(x(1) * (y(2) - y(3)) + x(2) * (y(3) - y(1)) + x(3) * (y(1) - y(2)));

H = (1/(2*A)).*[y(2) - y(3) , 0 , y(3) - y(1) , 0 , y(1) - y(2) , 0;
                x(3) - x(2) , 0 , x(1) - x(3) , 0 , x(2) - x(1) , 0;
                0 , y(2) - y(3) , 0 , y(3) - y(1) , 0 , y(1) - y(2);
                0 , x(3) - x(2) , 0 , x(1) - x(3) , 0 , x(2) - x(1)];
            
B0 = [y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0;
      0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1);
      x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),x(2)-x(1)];

ef = H * ed;
a = ef;
Au = [a(1),0,a(3),0;
      0,a(2),0,a(4);
      a(1),a(2),a(3),a(4)];

Bu = Au * H;
ee = (B0 + 0.5.*Bu)*ed;
end

