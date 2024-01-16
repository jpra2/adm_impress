function [ O, P, T, Qo ] = OPT_Interp_LPEW(ni)
global esurn1 esurn2 coord nsurn1 nsurn2 elem
%Retorna os vetores O, P, T e Qo.
% Lembrando que estes esurn1, nsurn1 j� estan ordenados em sentido
% anti-horario, sequencialmente. 

%Pr�-aloca��o dos vetores.%

P=zeros(nsurn2(ni+1)-nsurn2(ni),3); % vetor de pontos na vizinhan�a do n� "ni".
T=zeros(nsurn2(ni+1)-nsurn2(ni),3); % vetor de pontos dinamicos na vizinhan�a do n� "ni".
O=zeros(esurn2(ni+1)-esurn2(ni),3); % vetor de baricentro na vizinhan�a do n� "ni".
Qo=coord(ni,:);                     % coordenada do n� "ni".

%Constru��o dos vetores P, dos n�s vizinhos ao n� "ni", e T, dos pontos%
%m�dios das fases que concorrem no n� "ni".                            %

for i=1:size(P,1),
    P(i,:)=coord(nsurn1(nsurn2(ni)+i),:);
    T(i,:)=(P(i,:)+Qo)/2;
end

%Constru��o do vetor O, dos centr�ides (pontos de coloca��o) dos elementos%
%que concorrem no n� ni.                                                  %

for i=1:size(O,1),
    %Verifica se o elemento � um quadril�tero ou um tri�ngulo.
    if elem(esurn1(esurn2(ni)+i),4)==0 % lenbrando que o quarta columna
        b=3;                  
    else
        b=4;  % da matriz de elementos � para quadrilateros
    end
    %Carrega adequadamente o vetor O (braicentro de cada elemento)
    for j=1:b
        O(i,1)=O(i,1)+(coord(elem(esurn1(esurn2(ni)+i),j),1)/b);
        O(i,2)=O(i,2)+(coord(elem(esurn1(esurn2(ni)+i),j),2)/b);
        O(i,3)=O(i,3)+(coord(elem(esurn1(esurn2(ni)+i),j),3)/b);
    end
end

end
