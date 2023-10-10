//--------------------------------------------------------------------
//---------------- UNIVERSIDADE FEDERAL DE PERNAMBUCO ----------------
//---------------- CENTRO DE TECNOLOGIA E GEOCIENCIAS ----------------
//---------- PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL -----------
//--------------------------------------------------------------------

//Work developed by: Marcio Souza and Luiz E. Queiroz
//Adviser Professors: Paulo Lyra & Darlan Carvalho
//Create date: 2022/3/30;	hour: 8:37h

//--------------------------------------------------------------------
//This file has CAD parameters. It is related to building of domain

cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 0.5, 0, 1};
Point(4) = {1, 1, 0, 1};
Point(5) = {0, 1, 0, 1};
Point(6) = {0, 0.5, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {6, 3};
Line Loop(9) = {5, 7, 3, 4};
Plane Surface(9) = {9};
Line Loop(11) = {7, -2, -1, -6};
Plane Surface(11) = {11};
Physical Point(101) = {1, 5, 6};
Physical Point(102) = {2, 3, 4};
Physical Line(201) = {1, 4};
Physical Line(101) = {5, 6};
Physical Line(102) = {2, 3};
Physical Line(17) = {7};
Physical Surface(1) = {9, 11};


Transfinite Line {1,7} = 1 Using Progression 1.000000;
Transfinite Line {4,7} = 1 Using Progression 1.000000;
Transfinite Line {2,6} = 1 Using Progression 1.000000;
Transfinite Line {3,5} = 1 Using Progression 1.000000;
Transfinite Surface {9} = {3,4,5,6};

Recombine Surface {9};
