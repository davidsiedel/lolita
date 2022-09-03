//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.25, 0, 0, 1.0};
//+
Point(3) = {0.75, 0, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Point(5) = {1, 0.25, 0, 1.0};
//+
Point(6) = {1, 0.75, 0, 1.0};
//+
Point(7) = {1, 1, 0, 1.0};
//+
Point(8) = {0.75, 1, 0, 1.0};
//+
Point(9) = {0.25, 1, 0, 1.0};
//+
Point(10) = {0, 1, 0, 1.0};
//+
Point(11) = {0, 0.75, 0, 1.0};
//+
Point(12) = {0, 0.25, 0, 1.0};
//+
Point(13) = {0.25, 0.25, 0, 1.0};
//+
Point(14) = {0.75, 0.25, 0, 1.0};
//+
Point(15) = {0.75, 0.75, 0, 1.0};
//+
Point(16) = {0.25, 0.75, 0, 1.0};
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
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 1};
//+
Line(13) = {12, 13};
//+
Line(14) = {13, 2};
//+
Line(15) = {13, 14};
//+
Line(16) = {14, 3};
//+
Line(17) = {14, 5};
//+
Line(18) = {14, 15};
//+
Line(19) = {15, 6};
//+
Line(20) = {15, 8};
//+
Line(21) = {15, 16};
//+
Line(22) = {16, 9};
//+
Line(23) = {16, 11};
//+
Line(24) = {16, 13};
//+
Curve Loop(1) = {1, -14, -13, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {14, 2, -16, -15};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {16, 3, 4, -17};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {17, 5, -19, -18};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {19, 6, 7, -20};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {21, 22, -8, -20};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {23, -10, -9, -22};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {13, -24, 23, 11};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {15, 18, 21, 24};
//+
Plane Surface(9) = {9};
//+
Physical Point("NODE", 25) = {1};
//+
Physical Curve("BOTTOM", 26) = {1, 2, 3};
//+
Physical Curve("RIGHT", 27) = {4, 5, 6};
//+
Physical Curve("TOP", 28) = {7, 8, 9};
//+
Physical Curve("LEFT", 29) = {10, 11, 12};
//+
Physical Surface("SQUARE", 30) = {7, 6, 5, 4, 3, 2, 1, 8, 9};
//+
Transfinite Curve {23, 22, 9, 10, 20, 7, 6, 19, 17, 16, 3, 4, 12, 13, 14, 1} = 10 Using Progression 1;
//+
Transfinite Curve {11, 24, 21, 8, 18, 5, 15, 2} = 20 Using Progression 1;
