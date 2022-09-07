//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 0.5, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Point(6) = {0, 0.5, 0, 1.0};
//+
Point(7) = {0.3, 0.5, 0, 1.0};
//+
Point(8) = {0.7, 0.5, 0, 1.0};
//+
Point(9) = {0.5, 0.5, 0, 1.0};
//+
// Circle(1) = {7, 9, 8};
// //+
// Circle(2) = {8, 9, 7};
// //+
// Line(3) = {1, 2};
// //+
// Line(4) = {2, 3};
// //+
// Line(5) = {3, 4};
// //+
// Line(6) = {4, 5};
// //+
// Line(7) = {5, 6};
// //+
// Line(8) = {6, 1};
// //+
// Line(9) = {6, 7};
// //+
// Line(10) = {8, 3};
// //+
// Curve Loop(1) = {6, 7, 9, -2, 10, 5};
// //+
// Plane Surface(1) = {1};
// //+
// Curve Loop(2) = {8, 3, 4, -10, -1, -9};
// //+
// Plane Surface(2) = {2};//+
// Physical Curve("CIRCLE", 11) = {2, 1};
// //+
// Physical Curve("TOP", 12) = {6};
// //+
// Physical Curve("BOTTOM", 13) = {3};
// //+
// Physical Curve("LEFT", 14) = {7, 8};
// //+
// Physical Curve("RIGHT", 15) = {5, 4};
// //+ TOP
// Transfinite Curve {6} = 200 Using Progression 1;
// //+ TOP CIRCLE
// Transfinite Curve {2} = 120 Using Progression 1;
// //+ MIDDLE
// Transfinite Curve {9, 10} = 50 Using Progression 1;
// //+ TOP SIDES
// Transfinite Curve {7, 5} = 100 Using Progression 1;
// //+ BOTTOM CIRCLE
// Transfinite Curve {1, 1} = 50 Using Progression 1;
// //+ BOTTOM SIDES
// Transfinite Curve {8, -4} = 20 Using Progression 1.2;
// //+ BOTTOM
// Transfinite Curve {3} = 10 Using Progression 1;
// //+
// // Recombine Surface {1, 2};
//+
Point(10) = {0.5, 1, 0, 1.0};
//+
Point(11) = {0.5, 0.7, 0, 1.0};
//+
Point(12) = {0.5, 0.3, 0, 1.0};
//+
Point(13) = {0.5, 0, 0, 1.0};
//+
Line(1) = {1, 13};
//+
Line(2) = {13, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 10};
//+
Line(6) = {10, 5};
//+
Line(7) = {5, 6};
//+
Line(8) = {6, 1};
//+
Line(9) = {6, 7};
//+
Line(10) = {12, 13};
//+
Line(11) = {8, 3};
//+
Line(12) = {11, 10};
//+
Circle(13) = {7, 9, 11};
//+
Circle(14) = {11, 9, 8};
//+
Circle(15) = {8, 9, 12};
//+
Circle(16) = {12, 9, 7};
//+
Physical Curve("RIGHT", 17) = {4, 3};
//+
Physical Curve("BOTTOM", 18) = {1, 2};
//+
Physical Curve("LEFT", 19) = {8, 7};
//+
Physical Curve("TOP", 20) = {6, 5};
//+
Physical Curve("CIRCLE", 21) = {13, 14, 15, 16};
//+
Curve Loop(1) = {7, 9, 13, 12, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, -12, 14, 11, 4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, -11, 15, 10, 2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1, -10, 16, -9, 8};
//+
Plane Surface(4) = {4};
//+
Physical Surface("PLATE", 22) = {1, 2, 3, 4};
//+ TOP
Transfinite Curve {6, -5} = 40 Using Progression 1.02;
//+ TOP SIDES
Transfinite Curve {-7, 4} = 60 Using Progression 1.02;
//+ TOP VERT
Transfinite Curve {12} = 60 Using Progression 1.02;
//+ TOP CIRCLE
Transfinite Curve {13, 14} = 100 Using Progression 1;
//+ HORIZ
Transfinite Curve {9, 11} = 80 Using Progression 1;
//+ BOTTOM
Transfinite Curve {1, 2} = 10 Using Progression 1;
//+ BOTTOM SIDES
Transfinite Curve {8, -3} = 20 Using Progression 1.1;
//+ BOTTOM VERT
Transfinite Curve {10} = 10 Using Progression 1.1;
//+ BOTTOM CIRCLE
Transfinite Curve {-16, 15} = 20 Using Progression 1.1;
//+ TOP
// Transfinite Curve {6, -5} = 2 Using Progression 1.02;
// //+ TOP SIDES
// Transfinite Curve {-7, 4} = 2 Using Progression 1.02;
// //+ TOP VERT
// Transfinite Curve {12} = 2 Using Progression 1.02;
// //+ TOP CIRCLE
// Transfinite Curve {13, 14} = 2 Using Progression 1;
// //+ HORIZ
// Transfinite Curve {9, 11} = 2 Using Progression 1;
// //+ BOTTOM
// Transfinite Curve {1, 2} = 2 Using Progression 1;
// //+ BOTTOM SIDES
// Transfinite Curve {8, -3} = 2 Using Progression 1.1;
// //+ BOTTOM VERT
// Transfinite Curve {10} = 2 Using Progression 1.1;
// //+ BOTTOM CIRCLE
// Transfinite Curve {-16, 15} = 2 Using Progression 1.1;
//+
Physical Surface("SQUARE", 23) = {1, 2, 3, 4};
