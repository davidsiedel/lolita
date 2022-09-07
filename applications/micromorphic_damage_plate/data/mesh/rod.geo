//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 0.499, 0, 1.0};
//+
Point(3) = {0, 0.501, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0.1, 1, 0, 1.0};
//+
Point(6) = {0.1, 0.501, 0, 1.0};
//+
Point(7) = {0.1, 0.499, 0, 1.0};
//+
Point(8) = {0.1, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 1};
//+
Line(9) = {7, 2};
//+
Line(10) = {6, 3};
//+
Curve Loop(1) = {8, 1, -9, 7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 9, 2, -10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, 10, 3, 4};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {7, -1, -5, 3} = 50 Using Progression 1.05;
//+
Transfinite Curve {6, 2} = 2 Using Progression 1;
//+
Transfinite Curve {8, 4, 10, 9} = 6 Using Progression 1;
//+
Recombine Surface {1, 2, 3};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Physical Curve("BOTTOM", 11) = {8};
//+
Physical Curve("RIGHT", 12) = {7, 6, 5};
//+
Physical Curve("TOP", 13) = {4};
//+
Physical Curve("LEFT", 14) = {3, 2, 1};
//+
Physical Surface("SANE", 15) = {1, 3};
//+
Physical Surface("DEFECT", 16) = {2};
//+
Physical Surface("ROD", 17) = {3, 2, 1};
