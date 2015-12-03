Complete!

Modify and calfem_root.m so that the script can access the solveq.m script.
Or simply put the script in the same script as nr_solver_*.m

The plan3* methods were all written and verified, run test_plan_methods to
test.

a) The data.m script was written to show the mesh and pints where boundary
   conditions are set, with x (green) and y (red).
b) The nrchap6.m algorithm solves the system with a total L formulation,
   and gives good results.
c) I don't know which force is reffered to when considering the
   force-displacement curve, as we dont have an externa? force to work,
   with but it doesn't feel too important.