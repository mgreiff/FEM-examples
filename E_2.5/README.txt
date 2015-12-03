Modify and calfem_root.m so that the script can access the solveq.m script.
Or simply put the script in the same script as nr_solver_*.m

a) The three scripts can be tested by running test_bar_methods.m
b) The newtC2E5c.m was implemented and checked with lab assistants.
c) Run newtC2E5c.m to show the displacement as a function if the load.
d) Run newtC2E5d.m. The spring has been added to in the x-direction, but
   does not give good results - they are very similar to that of c), and I
   suspect we should get a continuous curve, as in 1.5.
e) Run newtC2E5e.m. Gives good results when the force is applied to the
   z-direction, but not when applied to the y-direction, as specified in
   the tsk - check with lab assistants if it should behave this way.
f) Run newtC2E5f.m. Will not work whn c = 0. Why?
g) Not yet written.