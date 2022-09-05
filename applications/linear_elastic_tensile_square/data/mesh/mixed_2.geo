//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 1.0};
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
Line(7) = {2, 5};
//+
Curve Loop(1) = {6, 1, 7, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, -4, -3, -2};
//+
Plane Surface(2) = {2};
//+
Physical Point("NODE", 8) = {1};
//+
Physical Curve("BOTTOM", 9) = {1, 2};
//+
Physical Curve("RIGHT", 10) = {3};
//+
Physical Curve("TOP", 11) = {4, 5};
//+
Physical Curve("LEFT", 12) = {6};
//+
Physical Surface("SQUARE", 13) = {1, 2};
//+
Transfinite Curve {6, 3, 7} = 50 Using Progression 1.01;
//+
Transfinite Curve {5, 4, 1, 2} = 25 Using Progression 1;
//+
Recombine Surface {2};
