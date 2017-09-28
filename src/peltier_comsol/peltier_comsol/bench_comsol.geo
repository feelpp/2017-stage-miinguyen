


Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthExtendFromBoundary=1;
Mesh.CharacteristicLengthFromPoints=1;
Mesh.ElementOrder=1;
Mesh.SecondOrderIncomplete = 0;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;
Mesh.RecombinationAlgorithm=0;
//Mesh.OptimizeNetgen=1;
// refine level 0
// partitioning data
/*
Mesh.Partitioner=1;
Mesh.NbPartitions=6;
Mesh.MshFilePartitioned=0;
*/

scaling=-3;
h=0.1*10^scaling;
Mesh.CharacteristicLengthMin=0;
Mesh.CharacteristicLengthMax=h;

ptBaseGeo0X=0;
ptBaseGeo0Y=0;
ptBaseGeo0Z=0;
lenghtGeo0X=1*10^scaling;
lenghtGeo0Y=1*10^scaling;
lenghtGeo0Z=0.1*10^scaling;
ptPairGeo0X=ptBaseGeo0X+lenghtGeo0X;
ptPairGeo0Y=ptBaseGeo0Y+lenghtGeo0Y;
ptPairGeo0Z=ptBaseGeo0Z+lenghtGeo0Z;

ptBaseGeo1X=ptBaseGeo0X;
ptBaseGeo1Y=ptBaseGeo0Y;
ptBaseGeo1Z=ptPairGeo0Z;
lenghtGeo1X=lenghtGeo0X;
lenghtGeo1Y=lenghtGeo0Y;
lenghtGeo1Z=5.8*10^scaling;
ptPairGeo1X=ptBaseGeo1X+lenghtGeo1X;
ptPairGeo1Y=ptBaseGeo1Y+lenghtGeo1Y;
ptPairGeo1Z=ptBaseGeo1Z+lenghtGeo1Z;

ptBaseGeo2X=ptBaseGeo0X;
ptBaseGeo2Y=ptBaseGeo0Y;
ptBaseGeo2Z=ptPairGeo1Z;
lenghtGeo2X=lenghtGeo0X;
lenghtGeo2Y=lenghtGeo0Y;
lenghtGeo2Z=0.1*10^scaling;
ptPairGeo2X=ptBaseGeo2X+lenghtGeo2X;
ptPairGeo2Y=ptBaseGeo2Y+lenghtGeo2Y;
ptPairGeo2Z=ptBaseGeo2Z+lenghtGeo2Z;



Point(1) = {ptBaseGeo1X,ptBaseGeo1Y,ptBaseGeo1Z, h};
Point(2) = {ptPairGeo1X,ptBaseGeo1Y,ptBaseGeo1Z, h};
Point(3) = {ptPairGeo1X,ptPairGeo1Y,ptBaseGeo1Z, h};
Point(4) = {ptBaseGeo1X,ptPairGeo1Y,ptBaseGeo1Z, h};
Point(5) = {ptBaseGeo1X,ptBaseGeo1Y,ptPairGeo1Z, h};
Point(6) = {ptPairGeo1X,ptBaseGeo1Y,ptPairGeo1Z, h};
Point(7) = {ptPairGeo1X,ptPairGeo1Y,ptPairGeo1Z, h};
Point(8) = {ptBaseGeo1X,ptPairGeo1Y,ptPairGeo1Z, h};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};
Line Loop(1) = {1,2,3,4};
Line Loop(2) = {-5,-8,-7,-6};
Line Loop(3) = {-1,9,5,-10};
Line Loop(4) = {10,6,-11,-2};
Line Loop(5) = {11,7,-12,-3};
Line Loop(6) = {8,-9,-4,12};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};


Point(9) = {ptBaseGeo0X,ptBaseGeo0Y,ptBaseGeo0Z, h};
Point(10) = {ptPairGeo0X,ptBaseGeo0Y,ptBaseGeo0Z, h};
Point(11) = {ptPairGeo0X,ptPairGeo0Y,ptBaseGeo0Z, h};
Point(12) = {ptBaseGeo0X,ptPairGeo0Y,ptBaseGeo0Z, h};
Line(13) = {9, 1};
Line(14) = {10, 2};
Line(15) = {9, 10};
Line(16) = {12, 4};
Line(17) = {12, 9};
Line(18) = {3, 11};
Line(19) = {11, 12};
Line(20) = {10, 11};

Line Loop(7) = {1, -14, -15, 13};
Line Loop(8) = {2, 18, -20, 14};
Line Loop(9) = {3, -16, -19, -18};
Line Loop(10) = {4, -13, -17, 16};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};

Point(13) = {ptBaseGeo2X,ptBaseGeo2Y,ptPairGeo2Z, h};
Point(14) = {ptPairGeo2X,ptBaseGeo2Y,ptPairGeo2Z, h};
Point(15) = {ptPairGeo2X,ptPairGeo2Y,ptPairGeo2Z, h};
Point(16) = {ptBaseGeo2X,ptPairGeo2Y,ptPairGeo2Z, h};
Line(21) = {5, 13};
Line(23) = {8, 16};
Line(24) = {7, 15};
Line(25) = {6, 14};
Line(26) = {13, 16};
Line(27) = {13, 14};
Line(28) = {14, 15};
Line(29) = {15, 16};

Line Loop(11) = {7, 23, -29, -24};
Plane Surface(11) = {11};
Line Loop(12) = {6, 24, -28, -25};
Plane Surface(12) = {12};
Line Loop(13) = {5, 25, -27, -21};
Plane Surface(13) = {13};
Line Loop(14) = {8, 21, 26, -23};
Plane Surface(14) = {14};
Line Loop(15) = {29, -26, 27, 28};
Plane Surface(15) = {15};
Line Loop(16) = {17, 15, 20, 19};
Plane Surface(16) = {16};

//materiaux centrale
Surface Loop(1) = {4, 3, 6, 5, 1, 2};
Volume(1) = {1};
//electrodes
Surface Loop(2) = {16, 10, 7, 8, 9, 1};
Volume(2) = {2};
Surface Loop(3) = {2, 15, 11, 14, 13, 12};
Volume(3) = {3};


//Physical Surface("mat0-Adiabatic") = {3,4,5,6};
Physical Surface("Ground") = {16};
Physical Surface("Intensity") = {15};
Physical Surface("Adiabatic") = {3,4,5,6,7,8,9,10,11,12,13,14};
Physical Surface("Internal-Boundaries") = {1,2};

Physical Volume("Material0") = {1};
Physical Volume("Material1") = {2,3};
//+
Physical Volume("Electrode1") = {3};
//+
Physical Volume("Electrode2") = {2};
