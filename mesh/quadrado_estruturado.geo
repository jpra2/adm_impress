// ==========================================
// 1. Definição da Geometria (Retângulo 6x6)
// ==========================================


Point(1) = {0, 0, 0, 1.0};
Point(2) = {10, 0, 0, 1.0};
Point(3) = {10, 10, 0, 1.0};
Point(4) = {0, 10, 0, 1.0};

Line(1) = {1, 2}; // Base
Line(2) = {2, 3}; // Direita
Line(3) = {3, 4}; // Topo
Line(4) = {4, 1}; // Esquerda

Line Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};

// ==========================================
// 2. Definição da Malha Estruturada
// ==========================================
// Regra de ouro: Para ter N elementos, você precisa de N+1 nós.
// Como queremos 6 quadrados de tamanho 1, precisamos de 7 pontos por linha.

Transfinite Line {1, 3} = 11; // 7 nós nas linhas horizontais (X)
Transfinite Line {2, 4} = 11; // 7 nós nas linhas verticais (Y)

// Aplica a interpolação estruturada na superfície
Transfinite Surface {1};

// Força o Gmsh a gerar quadriláteros (quads) em vez de triângulos
Recombine Surface {1};

// ==========================================
// 3. Grupos Físicos (Opcional, mas recomendado)
// ==========================================
// Ajuda a identificar as fronteiras e o domínio no seu solver (FEM/FVM)
// Physical Line("inferior") = {1};
// Physical Line("superior") = {3};
// Physical Line("esquerda") = {4};
// Physical Line("direita") = {2};
// Physical Surface("dominio") = {1};