
// 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
Mesh.Algorithm = 2;
Mesh.CharacteristicLengthFactor = 6;

k=1;

//divD = 30;
//divS = 30;
//divL = 10;
//nn = 1;
//diva =  18/nn;
//divb  = 15/nn;
//divc = 10/nn;
//divd = 10/nn;

// malha fina
nn = 0.193;

//malha grossa 2
//nn=1;

diva =  18/nn;
divb  = 15/nn;
divc = 10/nn;
divd = 10/nn;


cl__1 = 0.1;

x = 0.5*k;
y = 0.5*k;

a = 1*k;
b = 1*k;

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


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line {1} = 5 Using Progression 1;
Transfinite Line {2} = 5 Using Progression 1;
Transfinite Line {3} = 5 Using Progression 1;
Transfinite Line {4} = 5 Using Progression 1;


Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Recombine Surface {1};

Physical Line(101) = {4};
Physical Line(102) = {2};
Physical Line(201) = {1,3};
