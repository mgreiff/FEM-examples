clear all;close all;
% Checks that the plan3* methods have been written correctly by comparing to
% the Calfem computed values.

load('control6.1.mat')
[ compKe ] = plan3ge( ec' , t , D , ed' , es );
disp(compKe)
disp(Ke)
disp('---')
[ compef ] = plan3gf( ec', t, ed' , es );
disp(compef')
disp(ef)
disp(ef./compef')
disp('---')