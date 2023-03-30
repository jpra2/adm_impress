// Gmsh project created on Wed Mar 22 17:40:15 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 20, 0, 1.0};
//+
Point(4) = {0, 20, 0, 1.0};

//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Physical Line("Boundary") = {1,2,3,4};
Characteristic Length {1,2,3,4} = 1.0;
