// Points

Point(1) = {-2, -2, 0, 2};
Point(2) = {2, -2, 0, 2};
Point(3) = {2, 2, 0, 2};
Point(4) = {-2, 2, 0, 2};
Point(5) = {-10, -10, 0, 4};
Point(6) = {10, -10, 0, 4};
Point(7) = {10, 10, 0, 4};
Point(8) = {-10, 10, 0, 4};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Surfaces
Line Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};
Line Loop(2) = {7, 8, 5, 6};
Plane Surface(2) = {-1, -2};

// Meshing
Transfinite Line {3, 4, 1, 2} = 4 Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};



//+
Extrude {0, 0, 0.1} {
  Surface{1}; Surface{2}; Layers{5}; Recombine;
}
