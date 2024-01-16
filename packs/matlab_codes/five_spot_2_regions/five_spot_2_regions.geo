//--------------------------------------------------------------------
//---------------- UNIVERSIDADE FEDERAL DE PERNAMBUCO ----------------
//---------------- CENTRO DE TECNOLOGIA E GEOCIENCIAS ----------------
//---------- PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL -----------
//--------------------------------------------------------------------

//Work developed by: Marcio Souza and Luiz E. Queiroz
//Adviser Professors: Paulo Lyra & Darlan Carvalho
//Create date: 2022/4/19;	hour: 15:32h

//--------------------------------------------------------------------
//This file has CAD parameters. It is related to building of domain

//"cl1" corresponds to element size attributed in "Start.dat";
cl1 = 0.031300;

Point(1) = {0.000000, 0.000000, 0.000000, cl1};
Point(2) = {0.500000, 0.000000, 0.000000, cl1};
Point(3) = {1.000000, 0.000000, 0.000000, cl1};
Point(4) = {1.000000, 1.000000, 0.000000, cl1};
Point(5) = {0.500000, 1.000000, 0.000000, cl1};
Point(6) = {0.000000, 1.000000, 0.000000, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 5};
Line Loop(9) = {6, 1, 7, 5};
Plane Surface(9) = {9};
Line Loop(11) = {2, 3, 4, -7};
Plane Surface(11) = {11};
Physical Point(101) = {1, 2, 3, 4, 5, 6};
Physical Line(101) = {1, 2, 3, 4, 5, 6};
Physical Line(14) = {7};
Physical Surface(1) = {9, 11};

