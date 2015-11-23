Modify and calfem_root.m so that the script can access the solveq.m script.
Or simply put the script in the same script as nr_solver_*.m

The plan3* methods were all written and verified, run test_plan_methods to
test.

a) The data.m script was written to show the mesh and pints where boundary
   conditions are set, with x (green) and y (red).
b) The nrchap6.m algorithm solves the system with a total L formulation,
   but doesn't give good results. This, and the strange boundary conditions
   along x = 300 should be discussed during an exercise session.
c) This is done at the end of the script nrchap6.m.