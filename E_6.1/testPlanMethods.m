%clear all;close all;
% Checks that the plan3* methods have been written correctly by comparing to
% the Calfem computed values.

load('control6.1.mat')
[ ~ , compef ] = plan3gs( ec' , ed' );
disp(compef)
disp(ef)
disp('---')
[ compKe ] = plan3ge( ec' , t , D , ed' , es );
disp(compKe)
disp(Ke)
disp('---')