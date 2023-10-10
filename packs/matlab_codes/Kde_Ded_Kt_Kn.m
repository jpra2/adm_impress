function [Hesq, Kde, Kn, Kt, Ded] = Kde_Ded_Kt_Kn(kmap)

global bedge inedge coord elem centelem
% Retorna alguns parâmetros da expressão dos fluxos na face interna e no contorno.
% podemos verificar na pag. 5 e na pag. 6 do paper chinês

%Prealocação das matrizes.%

Hesq = zeros(size(bedge,1),1);
Kde = zeros(size(inedge,1),1);
Ded = zeros(size(inedge,1),1);
Kn = zeros(size(bedge,1),1);
Kt = zeros(size(bedge,1),1);
% determinar o flag do nó interior e fronteira de Neumann
K1 = zeros(3,3);
K2 = zeros(3,3);
K = zeros(3,3);

%Loop de arestas de contorno. ("bedge")
for ifacont = 1:size(bedge,1)
    %Determinação do baricentro a elemento à esquerda
    C1 = centelem(bedge(ifacont,3),:);
    %Determinação das alturas dos elementos à esquerda.
    ve1 = coord(bedge(ifacont,2),:) - coord(bedge(ifacont,1),:); % face
    ve2 = coord(bedge(ifacont,2),:) - C1; %Do centro esquerdo ao fim da face.
    
    ce = cross(ve1,ve2); % produto vetorial
    Hesq(ifacont) = norm(ce)/norm(ve1); % altura a relativo as faces do contorno
    
    
    %Essa é UMA maneira de construir os tensores
    K(1,1) = kmap(elem(bedge(ifacont,3),5),2);
    K(1,2) = kmap(elem(bedge(ifacont,3),5),3);
    K(2,1) = kmap(elem(bedge(ifacont,3),5),4);
    K(2,2) = kmap(elem(bedge(ifacont,3),5),5);
    
    %Cálculo das constantes tangenciais e normais
    
    Kn(ifacont) = (RotH(ve1)'*K*RotH(ve1))/norm(ve1)^2;
    Kt(ifacont) = (RotH(ve1)'*K*(ve1)')/norm(ve1)^2;
    
    %Adequação dos flags dos nós de Dirichlet
end  %End of FOR ("bedge")

%Loop de arestas internas. ("inedges")
for iface = 1:size(inedge,1),
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    C1 = centelem(inedge(iface,3),:); % baricentro do elemento a esquerda
    C2 = centelem(inedge(iface,4),:); % baricentro do elemento direito
    vcen = C2 - C1;
    vd1 = coord(inedge(iface,2),:) - coord(inedge(iface,1),:);
    
    %Determinação das alturas dos centróides dos elementos à direita e à%
    %esquerda.                                                          %
    
    vd2 = C2 - coord(inedge(iface,1),:);     %Do início da aresta até o
    %centro da célula da direita.
    cd = cross(vd1,vd2);
    H2 = norm(cd)/norm(vd1); % altura a direita
    
    ve2 = C1 - coord(inedge(iface,1),:);
    
    ce = cross(vd1,ve2);
    H1 = norm(ce)/norm(vd1); % altura a esquerda
    
    %Cálculo das constantes.%
    %A segunda entrada será tal que: 1=dir, 2=esq.
    
    %Essa é UMA maneira de construir os tensores.
    %Permeability on the Left
    
    K1(1,1) = kmap(elem(inedge(iface,3),5),2);
    K1(1,2) = kmap(elem(inedge(iface,3),5),3);
    K1(2,1) = kmap(elem(inedge(iface,3),5),4);
    K1(2,2) = kmap(elem(inedge(iface,3),5),5);
    
    %Permeability on the Right
    
    K2(1,1) = kmap(elem(inedge(iface,4),5),2);
    K2(1,2) = kmap(elem(inedge(iface,4),5),3);
    K2(2,1) = kmap(elem(inedge(iface,4),5),4);
    K2(2,2) = kmap(elem(inedge(iface,4),5),5);
    
    % calculo das constantes tangenciais e normais em cada face interna
    Kn1 = (RotH(vd1)'*K1*RotH(vd1))/norm(vd1)^2;
    Kt1 = (RotH(vd1)'*K1*(vd1)')/norm(vd1)^2;
    
    Kn2 = (RotH(vd1)'*K2*RotH(vd1))/norm(vd1)^2;
    Kt2 = (RotH(vd1)'*K2*(vd1)')/norm(vd1)^2;
    % calculo das constantes nas faces internas
    Kde(iface) = -norm(vd1)*((Kn1*Kn2))/(Kn1*H2 + Kn2*H1);
    % Ded: é uma constante que tem constantes geometricas + contantes
    % tangeciais
    Ded(iface) = (dot(vd1,vcen)/norm(vd1)^2) - ...
        (1/norm(vd1))*((Kt2/Kn2)*H1 + (Kt1/Kn1)*H2);
end  %End of FOR ("inedge")
end
function [RH]=RotH(vi)
%Função que retorna a rotação de um certo vetor em 90 graus, 
%anti-horariamente, na primeira coluna da matriz R, e horariamente, 
%na segunda. Restringe um vetor de 3 coordenadas a um de 2, considerando
%que a terceira coordenada é nula.
% vi2=zeros(2,1);
% if size(vi)~=[3 1]
%     vi=vi';
%     vi2(1)=vi(1);
%     vi2(2)=vi(2);
% end
RH=[0 1 0;-1 0 0;0 0 0]*vi';
end