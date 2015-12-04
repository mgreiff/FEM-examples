% Defines material properties
ex = [0, 10, 10;
      0, 10,  0];
  
ey = [0, 0,  10;
      0, 10, 10];

edof = [1 1 2 5 6 3 4;
        2 1 2 3 4 7 8];

E = 210*1e9;
v = 0.3;
t = 1; % OBS! arbitrarily set

nelm = length(edof(:,1));
ndof = max(max(edof));