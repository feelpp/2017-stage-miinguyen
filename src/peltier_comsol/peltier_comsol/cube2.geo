Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthExtendFromBoundary=1;
Mesh.CharacteristicLengthFromPoints=1;
Mesh.ElementOrder=1;
Mesh.SecondOrderIncomplete = 0;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;
//Mesh.OptimizeNetgen=1;
// partitioning data
Mesh.Partitioner=1;
Mesh.NbPartitions=4;
Mesh.MshFilePartitioned=0;
h = 0.1;
Mesh.RecombinationAlgorithm=0;
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
zmin = 0;
zmax = 1;
Point(1) = {xmin,ymin,zmin,h};
Point(2) = {xmax,ymin,zmin,h};
Point(3) = {xmax,ymax,zmin,h};
Point(4) = {xmin,ymax,zmin,h};
Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,1};
Line(4) = {1,2};
Point(5) = {xmax,ymin,zmax,h};
Point(6) = {xmin,ymin,zmax,h};
Line(5) = {2,5};
Line(6) = {5,6};
Line(7) = {6,1};
Point(7) = {xmin,ymax,zmax,h};
Point(8) = {xmax,ymax,zmax,h};
Line(8) = {6,7};
Line(9) = {5,8};
Line(10) = {8,7};
Line(11) = {4,7};
Line(12) = {3,8};



//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {2, 11, -10, -12};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {1, 12, -9, -5};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {4, 5, 6, 7};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {7, -3, 11, -8};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {6, 8, -10, -9};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 4, 3, 2, 5, 6};
//+
Volume(1) = {1};
//+
Physical Surface("Ground") = {2};
//+
Physical Surface("Intensity") = {4};
