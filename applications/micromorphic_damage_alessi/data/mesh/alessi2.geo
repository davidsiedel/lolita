//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 2, 0, 1.0};
//+
Point(4) = {0, 2, 0, 1.0};
//+
Point(5) = {0, 0.5, 0, 1.0};
//+
Point(6) = {1, 0.5, 0, 1.0};
//+
Point(7) = {1, 1.5, 0, 1.0};
//+
Point(8) = {0, 1.5, 0, 1.0};//+
Line(1) = {4, 3};
//+
Line(2) = {3, 7};
//+
Line(3) = {7, 8};
//+
Line(4) = {8, 4};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 2};
//+
Line(7) = {2, 1};
//+
Line(8) = {1, 5};
//+
Line(9) = {5, 8};
//+
Line(10) = {7, 6};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, -3, 10, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 5, 6, 7};
//+
Plane Surface(3) = {3};
//+ TIPS TOP
Transfinite Curve {1, 7} = 10 Using Progression 1;
//+ TIPS SIDES
Transfinite Curve {4, -2, -8, 6} = 15 Using Progression 1.2;
//+ TOPS
Transfinite Curve {3, 5} = 200 Using Progression 1;
//+ SIDES
Transfinite Curve {9, 10} = 200 Using Progression 1;
//+
// Transfinite Surface {2};
//+
Physical Curve("TOP", 11) = {1};
//+
Physical Curve("RIGHT", 12) = {2, 10, 6};
//+
Physical Curve("BOTTOM", 13) = {7};
//+
Physical Curve("LEFT", 14) = {8, 9, 4};
//+
Physical Surface("ROD", 15) = {1, 2, 3};
