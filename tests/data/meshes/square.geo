//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0.25, 0.25, 0, 1.0};
//+
Point(6) = {0.75, 0.25, 0, 1.0};
//+
Point(7) = {0.75, 0.75, 0, 1.0};
//+
Point(8) = {0.25, 0.75, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {8, 8};
//+
Line(8) = {7, 8};
//+
Line(9) = {8, 5};
//+
Curve Loop(1) = {8, 9, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(2) = {1, 2};
//+
Physical Point("NODE", 10) = {1};
//+
Physical Curve("TOP", 11) = {3};
//+
Physical Curve("BOTTOM", 12) = {1};
//+
Physical Curve("RIGHT", 13) = {2};
//+
Physical Curve("LEFT", 14) = {4};
//+
Transfinite Curve {3, 2, 1, 4} = 25 Using Progression 1;
//+
Transfinite Curve {8, 6, 5, 9} = 10 Using Progression 1;
//+
Recombine Surface {1};
//+
Transfinite Surface {1};
//+
Physical Surface("SQUARE", 15) = {1, 2};
