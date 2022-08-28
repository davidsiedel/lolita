//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("TOP", 5) = {3};
//+
Physical Curve("BOTTOM", 6) = {1};
//+
Physical Curve("LEFT", 7) = {4};
//+
Physical Curve("RIGHT", 8) = {2};
//+
Physical Surface("SQUARE", 9) = {1};
//+
Physical Point("NODE", 10) = {1};
//+
Transfinite Curve {3, 4, 1, 2} = 3 Using Progression 1;
//+
// Recombine Surface {1};
//+
// Transfinite Surface {1};
//+
// Recombine Surface {1};
//+
Transfinite Surface {1};
