SetFactory("OpenCASCADE");
h=0.1;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(7) = {6};
Physical Line("Intensity") = {1};
Physical Line("Ground") = {3};
Physical Surface(8) = {7};
//+
Physical Line("Material1") = {3, 1};
//+
Physical Surface("Material0") = {7};
