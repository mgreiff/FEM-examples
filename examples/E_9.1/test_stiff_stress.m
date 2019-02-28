load('controlDataC9E1A.mat')
disp('---')
Scalc = stresscal(ef,ep(1),ep(2));
disp(Scalc./es)
disp('---')
Dcalc = mstiff(ef,ep(1),ep(2));
disp(Dcalc./D)
disp('---')