h1 = 0.00002;
h2 = 0.015;
Point(1) = {0,0,0,h1};
Point(2) = {0,0.025,0,h1};
Point(3) = {0.0002,0,0,h1};
Point(4) = {0.0002,0.025,0,h1};
Point(5) = {0.5,0,0,h2};
Point(6) = {0.5,0.025,0,h2};
//+
Line(1) = {1, 2};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 2};
//+
Line(5) = {4, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {5, 3};
//+
Line Loop(1) = {-1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {-7, -6, -5, -3};
//+
Plane Surface(2) = {2};
//+
Physical Line(1) = {1};
//+
Physical Line(2) = {2, -7};
//+
Physical Line(3) = {6};
//+
Physical Line(4) = {-4, 5};
//+
Physical Line(5) = {3};
//+
Physical Surface(1) = {1};
//+
Physical Surface(2) = {2};
