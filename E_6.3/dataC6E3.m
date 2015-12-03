% Defines material properties
ex = [0, 10, 10;
      0, 10,  0];
  
ey = [0, 0,  10;
      0, 10, 10];

edof = [1 1 2 5 6 3 4;
        2 1 2 3 4 7 8];

E = 210*1e9;
v = 0.3;
t = 0.1; % OBS! arbitrarily set
D = (E/((1+v)*(1-2*v))) .* [1-v, v, 0;   v, 1-v, 0; 0, 0, (1-2*v)/2];

nelm = length(edof(:,1));
ndof = max(max(edof));