h = 0.1;
d = 0.4;
dd=d;
xmin = 0;
xmax=7.5+d;
ymin = 0;
ymax = 2;
// Box
Point(1) = {xmin-d, ymin-d, 0,h};
Point(2) = {xmax+d, ymin-d, 0,h};
Point(3) = {xmax+d, ymax+d, 0,h};
Point(4) = {xmin-d, ymax+d, 0,h};

// F
Point(30) = {0, 0, 0,h};
Point(31) = {0.2, 0, 0,h};
Point(32) = {0.2, 0.9, 0,h};
Point(33) = {0.8, 0.9, 0,h};
Point(34) = {0.8, 1.1, 0,h};
Point(35) = {0.2, 1.1, 0,h};
Point(36) = {0.2, 1.8, 0,h};
Point(37) = {1, 1.8, 0,h};
Point(38) = {1, 2, 0,h};
Point(39) = {0, 2, 0,h};

// E
Point(40) = {d+1, 0, 0,h};
Point(41) = {d+2, 0, 0,h};
Point(42) = {d+2, 0.2, 0,h};
Point(43) = {d+1.2, 0.2, 0,h};
Point(44) = {d+1.2, 0.9, 0,h};
Point(45) = {d+1.6, 0.9, 0,h};
Point(46) = {d+1.6, 1.1, 0,h};
Point(47) = {d+1.2, 1.1, 0,h};
Point(48) = {d+1.2, 1.8, 0,h};
Point(49) = {d+2,   1.8, 0,h};
Point(50) = {d+2,   2, 0,h};
Point(51) = {d+1,   2, 0,h};

// E
d=2*d;
Point(52) = {d+2, 0, 0,h};
Point(53) = {d+3, 0, 0,h};
Point(54) = {d+3, 0.2, 0,h};
Point(55) = {d+2.2, 0.2, 0,h};
Point(56) = {d+2.2, 0.9, 0,h};
Point(57) = {d+2.6, 0.9, 0,h};
Point(58) = {d+2.6, 1.1, 0,h};
Point(59) = {d+2.2, 1.1, 0,h};
Point(60) = {d+2.2, 1.8, 0,h};
Point(61) = {d+3,   1.8, 0,h};
Point(62) = {d+3,   2, 0,h};
Point(63) = {d+2,   2, 0,h};

// L
d=d+dd;
Point(64) = {d+3, 0, 0,h};
Point(65) = {d+4, 0, 0,h};
Point(66) = {d+4, 0.2, 0,h};
Point(67) = {d+3.2, 0.2, 0,h};
Point(74) = {d+3.2,   2, 0,h};
Point(75) = {d+3,   2, 0,h};

// +
d=d+dd;
Point(80) = {d+4,   0.6, 0,h};
Point(81) = {d+4.4, 0.6, 0,h};
Point(82) = {d+4.4, 0.2, 0,h};
Point(83) = {d+4.6, 0.2, 0,h};
Point(84) = {d+4.6, 0.6, 0,h};
Point(85) = {d+5.0, 0.6, 0,h};
Point(86) = {d+5.0, 0.8, 0,h};
Point(87) = {d+4.6, 0.8, 0,h};
Point(88) = {d+4.6, 1.2, 0,h};
Point(89) = {d+4.4, 1.2, 0,h};
Point(90) = {d+4.4, 0.8, 0,h};
Point(91) = {d+4  , 0.8, 0,h};

// +
d=d+dd;
Point(100) = {d+5,   0.6, 0,h};
Point(101) = {d+5.4, 0.6, 0,h};
Point(102) = {d+5.4, 0.2, 0,h};
Point(103) = {d+5.6, 0.2, 0,h};
Point(104) = {d+5.6, 0.6, 0,h};
Point(105) = {d+6.0, 0.6, 0,h};
Point(106) = {d+6.0, 0.8, 0,h};
Point(107) = {d+5.6, 0.8, 0,h};
Point(108) = {d+5.6, 1.2, 0,h};
Point(109) = {d+5.4, 1.2, 0,h};
Point(110) = {d+5.4, 0.8, 0,h};
Point(111) = {d+5  , 0.8, 0,h};

// box
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {30, 31};
Line(6) = {31, 32};
Line(7) = {32, 33};
Line(8) = {33, 34};
Line(9) = {34, 35};
Line(10) = {35, 36};
Line(11) = {36, 37};
Line(12) = {37, 38};
Line(13) = {38, 39};
Line(14) = {39, 30};
Line(15) = {40, 41};
Line(16) = {41, 42};
Line(17) = {42, 43};
Line(18) = {43, 44};
Line(19) = {44, 45};
Line(20) = {45, 46};
Line(21) = {46, 47};
Line(22) = {47, 48};
Line(23) = {48, 49};
Line(24) = {49, 50};
Line(25) = {50, 51};
Line(26) = {51, 40};
Line(27) = {63, 52};
Line(28) = {52, 53};
Line(29) = {53, 54};
Line(30) = {54, 55};
Line(31) = {55, 56};
Line(32) = {56, 57};
Line(33) = {57, 58};
Line(34) = {58, 59};
Line(35) = {59, 60};
Line(36) = {60, 61};
Line(37) = {61, 62};
Line(38) = {62, 63};
Line(39) = {64, 65};
Line(40) = {65, 66};
Line(41) = {66, 67};
Line(42) = {67, 74};
Line(43) = {74, 75};
Line(44) = {75, 64};
Line(45) = {82, 83};
Line(46) = {83, 84};
Line(47) = {84, 85};
Line(48) = {85, 86};
Line(49) = {86, 87};
Line(50) = {87, 88};
Line(51) = {88, 89};
Line(52) = {89, 90};
Line(53) = {90, 91};
Line(54) = {91, 80};
Line(55) = {80, 81};
Line(56) = {81, 82};
Line(57) = {102, 103};
Line(58) = {103, 104};
Line(59) = {104, 105};
Line(60) = {105, 106};
Line(61) = {106, 107};
Line(62) = {107, 108};
Line(63) = {108, 109};
Line(64) = {109, 110};
Line(65) = {110, 111};
Line(66) = {111, 100};
Line(67) = {100, 101};
Line(68) = {101, 102};
Line Loop(69) = {3, 4, 1, 2};
Line Loop(70) = {14, 5, 6, 7, 8, 9, 10, 11, 12, 13};
Line Loop(71) = {26, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
Line Loop(72) = {27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};
Line Loop(73) = {44, 39, 40, 41, 42, 43};
Line Loop(74) = {52, 53, 54, 55, 56, 45, 46, 47, 48, 49, 50, 51};
Line Loop(75) = {64, 65, 66, 67, 68, 57, 58, 59, 60, 61, 62, 63};
Plane Surface(76) = {69, 70, 71, 72, 73, 74, 75};

Physical Line("inlet") = {4};
Physical Line("wall") = {3, 1};
Physical Line("outlet") = {2};
Physical Line("letters") = {14, 13, 11, 12, 10, 9, 7, 8, 6, 5, 26, 24, 25, 23, 22, 21, 20, 19, 18, 17, 15, 16, 28, 30, 27, 32, 31, 34, 33, 35, 36, 37, 38, 43, 42, 44, 29, 39, 40, 41, 49, 47, 48, 50, 51, 52, 53, 54, 55, 56, 45, 46, 67, 66, 65, 63, 64, 62, 61, 60, 59, 58, 68, 57};
Physical Surface("feel") = {76};
