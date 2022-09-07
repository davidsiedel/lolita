//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, 0.71, 0, 1.0};
//+
Point(6) = {0.5, 0.7, 0, 1.0};
//+
Point(7) = {1, 0.71, 0, 1.0};
//+
Point(8) = {0.5, 0.3, 0, 1.0};
//+
Point(9) = {0.5, 0.5, 0, 1.0};
//+
Point(10) = {0.5, 1, 0, 1.0};
//+
Point(11) = {0.5, 0.71, 0, 1.0};
//+
Line(1) = {4, 10};
//+
Line(2) = {10, 3};
//+
Line(3) = {3, 7};
//+
Line(4) = {7, 2};
//+
Line(5) = {2, 1};
//+
Line(6) = {1, 5};
//+
Line(7) = {5, 4};
//+
Line(8) = {5, 11};
//+
Line(9) = {11, 7};
//+
Line(11) = {10, 11};
//+
Circle(12) = {6, 9, 8};
//+
Circle(13) = {8, 9, 6};
//+
Curve Loop(1) = {1, 11, -8, 7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, -9, -11};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 9, 4, 5, 6};
//+
Curve Loop(4) = {13, 12};
//+
Plane Surface(3) = {3, 4};
//+ TOP HORIZ
Transfinite Curve {1, 2, 9, 8} = 150 Using Progression 1;
//+ TOP VERTS
Transfinite Curve {7, -11, -3} = 50 Using Progression 1.1;
//+ BOTTOM VERTS
Transfinite Curve {-6, 4} = 25 Using Progression 1.15;
//+ CIRCLE
Transfinite Curve {12, -13} = 35 Using Progression 1.1;
//+ BOTTOM
Transfinite Curve {5} = 10 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
