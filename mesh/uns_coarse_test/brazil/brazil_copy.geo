
// 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
Mesh.Algorithm = 2;
Mesh.CharacteristicLengthFactor = 6;


//divD = 30;
//divS = 30;
//divL = 10;
//nn = 1;
//diva =  18/nn;
//divb  = 15/nn;
//divc = 10/nn;
//divd = 10/nn;


nn = 0.25;
diva =  18/nn;
divb  = 15/nn;
divc = 10/nn;
divd = 10/nn;


cl__1 = 0.1;

x = 0.5;
y = 0.5;

a = 1.5;
b = 1;

c = 1.2;
d = 0.75;

R = 0.2 ;
x = a/2;
y = b/2;


L = c/2;
l = d/2;
c = a/2;

//d = b/2 - 0.1;
//h = (b-a)/2;

//xi = x - (a/2);
//yi = y - (b/2);

Point(1) = {0, 0, 0, cl__1};
Point(2) = {a, 0, 0, cl__1};
Point(3) = {a, b, 0, cl__1};
Point(4) = {0, b, 0, cl__1};

Point(5) = {x, y - l, 0, cl__1};
Point(6) = {x - L, y , 0, cl__1};
Point(7) = {x , y  + l, 0, cl__1};
Point(8) = {x + L , y, 0, cl__1};


Point(9) = {x , y-R, 0, cl__1};

//Point(10) = {x+ R , y, 0, cl__1};

Point(10) = {x , y +R, 0, cl__1};
//Point(12) = {x -R , y, 0, cl__1};

Point(11) = {x , y, 0, cl__1};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


//Spline(5) = {5,6,7};
//Line(6) = {7,8};
//Spline(7) = {8,9,10};

//Line(8) = {10,11};
//Spline(9) = {11,12,13};
//Line(10) = {13,14};
//Spline(11) = {14,15,16};
//Line(12) = {16,5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

//BSpline(9) = {9, 10, 11,12,9};
//Spline(10) = {9, 11, 10};
//Spline(11) = {9, 11, 10};
//Spline(12) = {9, 11, 10};
Circle(9) = {9, 11, 10};
Circle(10) = {10, 11, 9};

Transfinite Line {1} = diva Using Progression 1;
Transfinite Line {2} = divb Using Progression 1;
Transfinite Line {3} = diva Using Progression 1;
Transfinite Line {4} = divb Using Progression 1;
Transfinite Line {5} = divc Using Progression 1;
Transfinite Line {6} = divc Using Progression 1;
Transfinite Line {7} = divc Using Progression 1;
Transfinite Line {8} = divc Using Progression 1;
Transfinite Curve {9} = divd Using Progression 1;
Transfinite Curve {10} = divd Using Progression 1;


//Transfinite Line {8} = divS Using Progression 1;
//Transfinite Line {12} = divL Using Progression 1;



Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Line Loop(3) = {9,10};

Plane Surface(1) = {1,2};
Plane Surface(2) = {2,3};

//Line {5} In Surface {2};
//Line {6} In Surface {2};
//Line {7} In Surface {2};
//Line {8} In Surface {2};
//Recombine Surface {1};
//Recombine Surface {2};

//Physical Point(101) = {1, 4};
//Physical Point(102) = {2, 3};
//Physical Point(103) = {9, 10};

//Physical Point(201) = {9,10};


//Physical Line(101) = {4};
//Physical Line(102) = {2};
//Physical Line(103) = {9,10};
//Physical Line(201) = {1,3};

Physical Line("Inflow") = {4};
Physical Line("Outflow") = {2};
Physical Line("InternalBoundary") = {9,10};
Physical Line("Walls") = {1,3};
Physical Line ("ExternalBoundary") = {1, 2, 3, 4};

//Physical Surface(1) = {1};
//Physical Surface(2) = {2};

Physical Surface("Green") = {1};
Physical Surface("Yellow") = {2};


