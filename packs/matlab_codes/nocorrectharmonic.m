function [y,weightDMP,raioaux]=nocorrectharmonic(kmap)

global inedge bedge coord elem centelem 
y=zeros(size(inedge,1)+size(bedge,1),3);
Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0;0 0 0];
raioaux=0;
for ifacont=1:size(bedge,1)
    
    y(ifacont,:)= 0.5*( coord(bedge(ifacont,1),:)+ coord(bedge(ifacont,2),:));
    
end
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    % calculo das projeções ortogonais sobre a face
    %[hrel]=projortogonal(rel,inedge(iface,2), inedge(iface,1));
    
    %[hlef]=projortogonal(lef,inedge(iface,1), inedge(iface,2));
    %Determinação das alturas dos centróides dos elementos
    
    %Do ponto do início da aresta até o centro da célula da direita
    vd2=centelem(rel,:)-coord(inedge(iface,1),:);
    cd=cross(vd1,vd2);
    hrel=norm(cd)/norm(vd1); % altura a direita
    
    %Do ponto do início da aresta até o centro da célula da direita
    ve2=centelem(lef,:)-coord(inedge(iface,1),:);
    ce=cross(vd1,ve2);
    hlef=norm(ce)/norm(vd1); % altura a esquerda
    
    % tensor do elemento a esquerda
    
    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);
    
    % tensor do elemento a direita
    
    Krel(1,1)=kmap(elem(rel,5),2);
    Krel(1,2)=kmap(elem(rel,5),3);
    Krel(2,1)=kmap(elem(rel,5),4);
    Krel(2,2)=kmap(elem(rel,5),5);
    
    % calculo das constantes normais em cada face interna
    Knlef=dot(R*vd1',Klef*(R*vd1')/norm(vd1)^2);
    
    Knrel=dot((R*(-vd1')),Krel*(R*(-vd1'))/norm(vd1)^2);
    % calculo dos pontos harmonicos 
    y(iface+size(bedge,1),:)=(hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)'+...
        hlef*hrel*(Klef'-Krel')*(R*(vd1/norm(vd1))'))/(hrel*Knlef+hlef*Knrel);
   
    
    % calculo dos pesos de interpolacao
    
    weightlef_1= (hrel*Knlef)/(hrel*Knlef+ hlef*Knrel); weightrel_1= 1-weightlef_1; % Eq. (17)
    weightlef_2 = 0;weightrel_2 = 0;
    weightlef = weightlef_1+weightlef_2;
    weightrel = weightrel_1+weightrel_2;
    weightDMP(iface,1)=weightlef;
    weightDMP(iface,2)=weightrel;
    weightDMP(iface,3)=lef;
    weightDMP(iface,4)=rel;
    %======================================================================%
end
end