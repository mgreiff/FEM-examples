clear all;
load('control2e5.mat')
compK = bar3ge( ec , ep , ed , es );
disp(compK)
disp(Ke)
disp('---')
[ ~ , N ] = bar3gs( ec , ep , ed );
disp(N)
disp(es)
disp('---')
[ force ] = bar3gf( ec , ed , es);
disp(force)
disp(ef)