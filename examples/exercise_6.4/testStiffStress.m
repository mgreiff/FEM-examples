load('data_e6.2.mat')
disp('---')
Scalc = stresscal(ef,ep(1),ep(2));
disp(Scalc./es)
disp('---')
Dcalc = mstiff(ef,ep(1),ep(2));
disp(Dcalc./D)
disp('---')