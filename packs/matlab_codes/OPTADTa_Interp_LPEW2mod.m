function [ O, P, T, Qo, areas, distancias, talfa] = OPTADTa_Interp_LPEW2mod(ni)
global esurn1 esurn2 coord nsurn1 nsurn2 elem
%Retorna os vetores O, P, T, Qo, areas, distancias e talfa.
%Lembrando que estes esurn1, nsurn1 já estão ordenados em sentido anti-horario, sequencialmente. 

%Pré-alocação dos vetores.%

P=zeros(nsurn2(ni+1)-nsurn2(ni),3); % vetor de pontos na vizinhança do nó "ni".
T=zeros(nsurn2(ni+1)-nsurn2(ni),3); % vetor de pontos dinamicos na vizinhança do nó "ni".
O=zeros(esurn2(ni+1)-esurn2(ni),3); % vetor de baricentro na vizinhança do nó "ni".
Qo=coord(ni,:);                     % coordenada do nó "ni".
distancias=zeros(esurn2(ni+1)-esurn2(ni),2); %vetor das distâncias ortogonais necessárias.
talfa=zeros(esurn2(ni+1)-esurn2(ni),2); %vetor dos vetores tangenciais.
areas=zeros(esurn2(ni+1)-esurn2(ni),2);%vetor das áreas.


%Construção dos vetores P, dos nós vizinhos ao nó "ni", e T, dos pontos%
%médios das fases que concorrem no nó "ni".                            %

for i=1:size(P,1),
    P(i,:)=coord(nsurn1(nsurn2(ni)+i),:);
    T(i,:)=(P(i,:)+Qo)/2;
    %Construção dos vetores tangenciais%
    epsilon=0.5; %definido a princípio%
    %talfa(i,1:2)=(tsigma,ttau)%;
    talfa(i,1)=(T(i,:)-Qo)/epsilon; %Esse é o tsigma%
    talfa(i,2)=(T(i+1,:)-Qo)/epsilon; %Esse é o ttau%
end


%Construção do vetor O, dos centróides (pontos de colocação) dos elementos%
%que concorrem no nó ni.                                                  %

for i=1:size(O,1),
    %Verifica se o elemento é um quadrilátero ou um triângulo.
    if elem(esurn1(esurn2(ni)+i),4)==0 % lenbrando que o quarta columna
        b=3;                  
    else
        b=4;  % da matriz de elementos é para quadrilateros
    end
    %Carrega adequadamente o vetor O (braicentro de cada elemento)
    for j=1:b
        O(i,1)=O(i,1)+(coord(elem(esurn1(esurn2(ni)+i),j),1)/b);
        O(i,2)=O(i,2)+(coord(elem(esurn1(esurn2(ni)+i),j),2)/b);
        O(i,3)=O(i,3)+(coord(elem(esurn1(esurn2(ni)+i),j),3)/b);
    end
end


%Obtenção da área da região triangular%
for i=1:size(O,1)
    %Calcula os vetores necessários%
    v1=T(i,:)-Qo;
    v2=T(i+1,:)-Qo;
    areas(i,1)=0.5*cross(v1,v2);
end


%Construção das distâncias ortogonais%
for i=1:size(O,1)

    %distancias(i,1:2)=(distanciasigma,distanciatau)%

%Preenchimentos das distancias tau%
    if(i==size(distancias,1))&&(size(P,1)==size(O,1))
       aux1=O(i,:)-Qo;
       aux2=P(1,:)-Qo;
       distancias(i,2)=norm(cross(aux1,aux2))/norm(aux2); %está é a d k,tau%
       %distancias(i-1,2)=norm(cross(aux1,aux2))/norm(aux2); %está é a d ksigma,tau%
       %distancias(i+1,2)=norm(cross(aux1,aux2))/norm(aux2); %está é a d ktau,tau%
    else
       aux1=O(i,:)-Qo;
       aux2=P(i+1,:)-Qo;
       distancias(i,2)=norm(cross(aux1,aux2))/norm(aux2); %está é a d k,tau%
       %distancias(i-1,2)=norm(cross(aux1,aux2))/norm(aux2); %está é a d ksigma,tau%
       %distancias(i+1,2)=norm(cross(aux1,aux2))/norm(aux2); %está é a d ktau,tau%
    end
  
%Preenchimentos das distancias sigma%
       aux1=O(i,:)-Qo;
       aux2=P(1,:)-Qo;
       distancias(i,1)=norm(cross(aux1,aux2))/norm(aux2); %está é a d k,sigma%
       %distancias(i-1,1)=norm(cross(aux1,aux2))/norm(aux2); %está é a d ksigma,sigma%
       %distancias(i+1,1)=norm(cross(aux1,aux2))/norm(aux2); %está é a d ksigma,sigma%

end
end
    
