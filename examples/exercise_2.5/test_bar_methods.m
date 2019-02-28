clear all;close all;
% Checks that the bar3* methods have been written correctly by comparing to
% the Calfem computed values. Note that we work with a transposed force
% vector.

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