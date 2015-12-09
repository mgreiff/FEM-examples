clear all;close all;
% Checks that the plan3* methods have been written correctly by comparing to
% the Calfem computed values.

load('controlPlan3Methods.mat')
[ compKe ] = plan3ge( ec' , t , D , ed' , es );
disp(compKe./Ke)
disp('---')
[ compef ] = plan3gf( ec', t, ed' , es );
disp(compef'./ef)
disp('---')