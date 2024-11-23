function [y,weightDMP,raioaux]=secondcorrectharmonic(kmap,N)

global inedge bedge coord elem centelem benchmark
y=zeros(size(inedge,1)+size(bedge,1),3);
Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0;0 0 0];
raioaux=0;
for ifacont=1:size(bedge,1)
    
    y(ifacont,:)= 0.5*( coord(bedge(ifacont,1),:)+ coord(bedge(ifacont,2),:));
    
end
contador=1;
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
    % calculo dos pontos armonicos
    y(iface+size(bedge,1),:)=(hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)'+...
        hlef*hrel*(Klef'-Krel')*(R*(vd1/norm(vd1))'))/(hrel*Knlef+hlef*Knrel);
    
    %=====================================================================%
    %     calculo dos pontos que saim fora da face interna nos pontos
    %     intermeadiarios da discontinuidade, sobre na discontinuidade
    %     horizontal.
    if strcmp(benchmark,'benchmar5_6') || strcmp(benchmark,'benchmar5_7')
        r1= find(N(inedge(iface,1),:)~=0 & N(inedge(iface,1),:)~=iface);
        a1=N(inedge(iface,1),r1);
        
        s1= find(N(inedge(iface,2),:)~=0 & N(inedge(iface,2),:)~=iface );
        b1=N(inedge(iface,2),s1);
        
        vetortotal= [a1 b1];
        
        for j=vetortotal
            
            if (j<size(inedge,1) || j==size(inedge,1))
                a1=inedge(j,1);
                a2=inedge(j,2);
                if (((coord(a1,1)< y(iface+size(bedge,1),1) || abs(coord(a1,1)- y(iface+size(bedge,1),1))<1e-10)  && (y(iface+size(bedge,1),1)< coord(a2,1) || abs(y(iface+size(bedge,1),1)- coord(a2,1))<1e-10)) &&...
                        ((coord(a1,2)< y(iface+size(bedge,1),2)|| abs(coord(a1,2)- y(iface+size(bedge,1),2))<1e-10) && (y(iface+size(bedge,1),2) < coord(a2,2) || abs(y(iface+size(bedge,1),2)-coord(a2,2))<1e-10))) ||...
                        (((coord(a2,1)< y(iface+size(bedge,1),1) || abs(coord(a2,1)- y(iface+size(bedge,1),1))<1e-10)  && (y(iface+size(bedge,1),1)< coord(a1,1) || abs(y(iface+size(bedge,1),1)- coord(a1,1))<1e-10)) &&...
                        ((coord(a2,2)< y(iface+size(bedge,1),2)|| abs(coord(a2,2)- y(iface+size(bedge,1),2))<1e-10) && (y(iface+size(bedge,1),2) < coord(a1,2) || abs(y(iface+size(bedge,1),2)-coord(a1,2))<1e-10)))
                    % levando ao ponto interior
                    y(iface+size(bedge,1),:)= (hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)')/(hrel*Knlef+hlef*Knrel);
                    
                    break
                end
            else
                a1=bedge(j-size(inedge,1),1);
                a2=bedge(j-size(inedge,1),2);
                if (((coord(a1,1)< y(iface+size(bedge,1),1) || abs(coord(a1,1)- y(iface+size(bedge,1),1))<1e-10)  && (y(iface+size(bedge,1),1)< coord(a2,1) || abs(y(iface+size(bedge,1),1)- coord(a2,1))<1e-10)) &&...
                        ((coord(a1,2)< y(iface+size(bedge,1),2)|| abs(coord(a1,2)- y(iface+size(bedge,1),2))<1e-10) && (y(iface+size(bedge,1),2) < coord(a2,2) || abs(y(iface+size(bedge,1),2)-coord(a2,2))<1e-10))) ||...
                        (((coord(a2,1)< y(iface+size(bedge,1),1) || abs(coord(a2,1)- y(iface+size(bedge,1),1))<1e-10)  && (y(iface+size(bedge,1),1)< coord(a1,1) || abs(y(iface+size(bedge,1),1)- coord(a1,1))<1e-10)) &&...
                        ((coord(a2,2)< y(iface+size(bedge,1),2)|| abs(coord(a2,2)- y(iface+size(bedge,1),2))<1e-10) && (y(iface+size(bedge,1),2) < coord(a1,2) || abs(y(iface+size(bedge,1),2)-coord(a1,2))<1e-10)))
                    
                    % levando ao ponto interior
                    y(iface+size(bedge,1),:)= (hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)')/(hrel*Knlef+hlef*Knrel);
                    
                    break
                end
            end
        end
    else
        % ative este calculo em malhas severamente distorcidas ja que
        a=0.5*(coord(inedge(iface,1),:)+coord(inedge(iface,2),:)); % ponto medio da face
        % este calculo eh usad por varios autores na literatura
        
        
        %=====================================================================%
        % comprimento da face ou area
        RR=0.5*norm(vd1);
        
        % R' is called as Raux
        Raux=3*RR;
        if norm(y(iface+size(bedge,1),:)- a)>Raux
            
            
            % if norm(a-coord(inedge(iface,1),:))<norm(a-y(iface+size(bedge,1),:))
            % correccao ao ponto medio
            y(iface+size(bedge,1),:)= a;
            contador=contador+1;
        end
        
    end
    
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
name=(contador/size(inedge,1))*100;
sprintf('>> Percentage of faces corrected: %s ',num2str(name))
end