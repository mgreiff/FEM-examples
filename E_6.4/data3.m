
gnew('geom')

gpoints(1,0,0,0)
gpoints(2,50,0,0)
gpoints(3,150,0,0)
gpoints(4,300,0,0)
gpoints(5,300,150,0)
gpoints(6,150,150,0)
gpoints(7,0,150,0)
gpoints(8,0,50,0)
gpoints(9,cos(pi/4)*50,sin(pi/4)*50,0)

glines(1,8,2,3)
glines(2,8,3,4)
glines(3,8,4,5)
glines(4,8,5,6)
glines(5,8,6,3)
glines(6,8,6,7)
glines(7,8,7,8)
glines(8,8,8,9,1)
glines(9,8,9,2,1)
glines(10,8,6,9)


gsurfs(1,2,3,4,5)
gsurfs(2,1,5,10,9)
gsurfs(3,10,6,7,8)
gmesh('surf',1,'qu4',2,1);
gmesh('surf',2,'qu4',2,1);
gmesh('surf',3,'qu4',2,1);

[edof,coord,dof]=getgeom;

[ex,ey,ez]=coordxtr(edof,coord,dof,4);

bc=gboundary('line',1,2,0);
bc=gboundary('line',2,2,0,bc);
bc=gboundary('line',7,1,0,bc);
bc=gboundary('line',3,1,10,bc);
