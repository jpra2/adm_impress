//Mesh.Algorithm = 2;

cl__1 = 1e+22;
cl__2 = 1e+22;

Lx = 100;
Ly = 100;

Point(1) = {0, 0, 0, cl__1};
Point(2) = {Lx, 0, 0, cl__1};
Point(3) = {Lx, Ly, 0, cl__1};
Point(4) = {0, Ly, 0, cl__1};



Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Curve {1, 2, 3, 4} = 8 + 1 Using Progression 1;
//Transfinite Curve {11, 9, 14, 12} = 10 + 1 Using Progression 1;
//Transfinite Surface {1};
Transfinite Surface {1};
//Transfinite Surface {3};
Recombine Surface{1};

//Line {5} In Surface {1};
//Line {6} In Surface {1};
//Line {7} In Surface {1};
//Line {8} In Surface {1};
//Line {9} In Surface {1};
//Line {10} In Surface {1};
//Line {11} In Surface {1};
//Line {12} In Surface {1};
//Physical Point(201) = {1, 2, 3, 4};
//Physical Curve(201) = {8};
//Physical Curve(202) = {4};
//Physical Curve(203) = {1, 2, 3, 5, 6, 7};
//Physical Surface(1) = {1};
//Physical Surface(2) = {2};
//Physical Surface(3) = {3};

//Mesh 2;
