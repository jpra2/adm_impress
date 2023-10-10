function [trasmi]=transmissibility(p,kmap,gamma,Fface,parameter,pinterp,valuemin,pointcolo)
global inedge coord elem bedge
Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0;0 0 0];
for ifacecont=1:size(bedge,1)
    norma= norm(coord(bedge(ifacecont,1),:)-coord(bedge(ifacecont,2),:));
        
    %% ej(iface) do elemento a esquerda
    
    trasmi(ifacecont,1)= norma*(parameter(1,2,ifacecont));
        
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    %Determinação das alturas dos centróides dos elementos à direita e à
    %esquerda.
    
    vd2=pointcolo(inedge(iface,4),:)-coord(inedge(iface,1),:);     %Do início da aresta até o
    %centro da célula da direita.
    cd=cross(vd1,vd2);
    hrel=norm(cd)/norm(vd1); % altura a direita
    %hrel= norm((dot(vd2,vd1)/(norm(vd1)^2))*vd1);
    ve2=pointcolo(inedge(iface,3),:)-coord(inedge(iface,1),:);
    
    ce=cross(vd1,ve2);
    hlef=norm(ce)/norm(vd1); % altura a esquerda
    %hlef= norm((dot(ve2,vd1)/(norm(vd1)^2))*vd1);
    %Cálculo das constantes.%
    %A segunda entrada será tal que: 1=dir, 2=esq.
    
    %Essa é UMA maneira de construir os tensores.
    
    Klef(1,1)=kmap(elem(inedge(iface,3),5),2);
    Klef(1,2)=kmap(elem(inedge(iface,3),5),3);
    Klef(2,1)=kmap(elem(inedge(iface,3),5),4);
    Klef(2,2)=kmap(elem(inedge(iface,3),5),5);
    
    % tensor a elemento a direita
    
    Krel(1,1)=kmap(elem(inedge(iface,4),5),2);
    Krel(1,2)=kmap(elem(inedge(iface,4),5),3);
    Krel(2,1)=kmap(elem(inedge(iface,4),5),4);
    Krel(2,2)=kmap(elem(inedge(iface,4),5),5);
    
    Knlef=((R*vd1')'*Klef*(R*vd1'))/norm(vd1)^2;
    
    Knrel=((R*(-vd1'))'*Krel*(R*(-vd1')))/norm(vd1)^2;
    %% calculo dos pesos
    % calculo do peso correspondente a face direita
    pesolef= (hrel*Knlef)/(hrel*Knlef+ hlef*Knrel);
    % calculo do peso correspondente a face esquerda
    pesorel= 1-pesolef;
    
    [mulef, murel,partialfluxlef2,partialfluxrel2]=calculomu(lef, rel,Fface,valuemin,p,pinterp,iface,gamma,parameter);
    
    % Calculo da transimissibilidade
    
    TT=(1-gamma)*norma*(mulef*parameter(1,2,iface+size(bedge,1))*pesorel + murel*parameter(2,2,iface+size(bedge,1))*pesolef);
    
    trasmi(iface+size(bedge,1),1)= TT;
end
end