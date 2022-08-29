//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {44, 0, 0, 1.0};
//+
Point(3) = {48, 60, 0, 1.0};
//+
Point(4) = {48, 44, 0, 1.0};
//+
Point(5) = {0, 44, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 5};
//+
Line(4) = {5, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("LEFT", 5) = {4};
//+
Physical Curve("RIGHT", 6) = {2};
//+
Physical Surface("SQUARE", 7) = {1};
//+
Transfinite Curve {3, 1} = 40 Using Progression 1;
//+
Transfinite Curve {4, 2} = 20 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
