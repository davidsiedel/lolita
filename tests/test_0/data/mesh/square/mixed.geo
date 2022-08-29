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
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {7, 8, 5, 6};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Point("NODE", 9) = {1};
//+
Physical Curve("TOP", 10) = {3};
//+
Physical Curve("BOTTOM", 11) = {1};
//+
Physical Curve("LEFT", 12) = {4};
//+
Physical Curve("RIGHT", 13) = {2};
//+
Physical Surface("SQUARE", 14) = {1, 2};
//+
// Transfinite Curve {3, 4, 1, 2, 6, 7, 8, 5} = 61 Using Progression 1;
//+
Transfinite Curve {3, 4, 1, 2} = 101 Using Progression 1;
//+
Transfinite Curve {6, 7, 8, 5} = 51 Using Progression 1;
//+
Recombine Surface {2};
