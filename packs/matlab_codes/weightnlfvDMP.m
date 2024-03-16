function [weightDMP]=weightnlfvDMP(kmap,y)
global inedge centelem elem coord bedge
Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0;0 0 0];
for iface=1:size(inedge,1);
    lef=inedge(iface,3);
    rel=inedge(iface,4);
vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
%% Determinação das alturas dos centróides dos elementos

vd2=centelem(inedge(iface,4),:)-coord(inedge(iface,1),:); %Do início da aresta até o
%centro da célula da direita.
cd=cross(vd1,vd2);
% altura a direita
hrel=norm(cd)/norm(vd1); 
ve2=centelem(inedge(iface,3),:)-coord(inedge(iface,1),:);

ce=cross(vd1,ve2);
% altura a esquerda
hlef=norm(ce)/norm(vd1); 
%Cálculo das constantes.%

% tensor do elemento esquerda

Klef(1,1)=kmap(elem(inedge(iface,3),5),2);
Klef(1,2)=kmap(elem(inedge(iface,3),5),3);
Klef(2,1)=kmap(elem(inedge(iface,3),5),4);
Klef(2,2)=kmap(elem(inedge(iface,3),5),5);

% tensor do elemento direita

Krel(1,1)=kmap(elem(inedge(iface,4),5),2);
Krel(1,2)=kmap(elem(inedge(iface,4),5),3);
Krel(2,1)=kmap(elem(inedge(iface,4),5),4);
Krel(2,2)=kmap(elem(inedge(iface,4),5),5);
a=0.5*(coord(inedge(iface,1),:)+coord(inedge(iface,2),:)); % ponto medio da face
% calculo das constantes tangenciais e normais em cada face interna
Knlef=((R*vd1')'*Klef*(R*vd1'))/norm(vd1)^2;

Knrel=((R*(-vd1'))'*Krel*(R*(-vd1')))/norm(vd1)^2;
% if (1/2)*norm(vd1)<norm(a-y(iface+size(bedge,1),:))
%    %Ponto harmônico corrigido:
%     weightlef_2= (Knlef/hlef)/(Knlef/hlef + Knrel/hrel); weightrel_2= 1-weightlef_2;
%     weightlef_1 = 0; weightrel_1 = 0;
% else 
    %Ponto harmônico inicial:
    weightlef_1= (hrel*Knlef)/(hrel*Knlef+ hlef*Knrel); weightrel_1= 1-weightlef_1; % Eq. (17)
     weightlef_2 = 0;weightrel_2 = 0;
% end
weightlef = weightlef_1+weightlef_2;
weightrel = weightrel_1+weightrel_2;

%% Calculo dos pesos
%weightlef= (hrel*Knlef)/(hrel*Knlef+ hlef*Knrel); weightrel= 1-weightlef; % Eq. (17)
weightDMP(iface,1)=weightlef;
weightDMP(iface,2)=weightrel;
weightDMP(iface,3)=lef;
weightDMP(iface,4)=rel;
end
end