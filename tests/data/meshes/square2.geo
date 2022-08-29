//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.4, 0, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Point(5) = {0.5, 1, 0, 1.0};
//+
Point(6) = {0, 1, 0, 1.0};
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
Line(6) = {6, 1};
//+
Line(7) = {5, 2};
//+
Curve Loop(1) = {5, 6, 1, -7};
//+
Curve Loop(2) = {4, 7, 2, 3};
//+
Plane Surface(1) = {1};
//+
Plane Surface(2) = {2};
//+
Physical Point("NODE", 8) = {1};
//+
Physical Curve("BOTTOM", 9) = {1, 2};
//+
Physical Curve("TOP", 10) = {5, 4};
//+
Physical Curve("LEFT", 11) = {6};
//+
Physical Curve("RIGHT", 12) = {3};
//+
Physical Surface("SQUARE", 13) = {1, 2};
//+
Recombine Surface {2};
//+
Transfinite Curve {5, 6, 1, 2, 7, 3, 4} = 2 Using Progression 1;
//+
Transfinite Surface {1};
