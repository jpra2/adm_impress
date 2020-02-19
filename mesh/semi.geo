k_pe_to_m = 1/3.28084;
k2 = 1.0;
Lx=900*k2;
Ly=450*k2;
Lz=90*k2;
hx = 20*k2;
hy = 10*k2;
hz = 2*k2;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {Lx, 0, 0, 1.0};
//+
Point(3) = {0, Ly, 0, 1.0};
//+
Point(4) = {Lx, Ly, 0, 1.0};
//+
Point(5) = {0, 0, Lz, 1.0};
//+
Point(6) = {Lx, 0, Lz, 1.0};
//+
Point(7) = {0, Ly, Lz, 1.0};
//+
Point(8) = {Lx, Ly, Lz, 1.0};
//+
Line(1) = {5, 6};
//+
Line(2) = {6, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 5};
//+
Line(5) = {5, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 4};
//+
Line(8) = {4, 3};
//+
Line(9) = {3, 7};
//+
Line(10) = {3, 1};
//+
Line(11) = {6, 8};
//+
Line(12) = {2, 4};
//+
Curve Loop(1) = {1, 11, -6, -5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, -9, 10, 4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, 7, 8, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {10, -3, 12, 8};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {3, 4, 1, 2};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {2, 12, -7, -11};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 5, 4, 2, 3, 6};

//+
Volume(1) = {1};

//+
Transfinite Curve {6, 1, 3, 8} = Lx/hx + 1 Using Progression 1;
Transfinite Curve {5, 10, 11, 12} = Ly/hy + 1 Using Progression 1;
Transfinite Curve {2, 4, 7, 9} = Lz/hz + 1 Using Progression 1;

//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {6};
//+
Transfinite Surface {5};
//+
Transfinite Volume {1};
//+
Recombine Surface {1, 3, 6, 5, 2, 4};
