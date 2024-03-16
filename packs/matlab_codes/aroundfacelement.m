function [ksii,ksij,aux11,aux12,auxy]=aroundfacelement(facelement,y,ielem,ve2,kielem,kmap,raioaux)
global bedge inedge centelem elem coord
% esta rutina reconstrui o estencil para evitar inconsistencia no calculo
% do fluxo.
Klef=zeros(3,3);
Krel=zeros(3,3);
s=facelement(ielem,1:kielem);
R=[0 1 0; -1 0 0;0 0 0];
for auxiface=s
    if auxiface<size(bedge,1) || auxiface==size(bedge,1)
        
        v1=bedge(auxiface,1);
        v2=bedge(auxiface,2);
        y(auxiface,:)=0.5*(coord(v1,:)+coord(v2,:));
    else
        v1=inedge(auxiface-size(bedge,1),1);
        v2=inedge(auxiface-size(bedge,1),2);
        lef=inedge(auxiface-size(bedge,1),3);
        rel=inedge(auxiface-size(bedge,1),4);
         
         vd1=coord(v2,:)-coord(v1,:);
        %Determinação das alturas dos centróides dos elementos
        
        %Do ponto do início da aresta até o centro da célula da direita
         vd2=centelem(rel,:)-coord(v1,:);
         cd=cross(vd1,vd2);
         hrel=norm(cd)/norm(vd1); % altura a direita
%         
%         %Do ponto do início da aresta até o centro da célula da direita
         ve21=centelem(lef,:)-coord(v1,:);
         ce=cross(vd1,ve21);
         hlef=norm(ce)/norm(vd1); % altura a esquerda
         
%         % tensor do elemento a esquerda
         
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
         
        % y(auxiface,:)=(hrel*Knlef*centelem(lef,:)'+ ...
        % hlef*Knrel*centelem(rel,:)')/(hrel*Knlef+hlef*Knrel);  % sempre usando 
        %y(auxiface,:)=((Knlef/hlef)*coord(v1,:)'+(Knrel/hrel)*coord(v2,:)')/((Knlef/hlef)+(Knrel/hrel));
         y(auxiface,:)=0.5*(coord(v1,:)+coord(v2,:));
    end
    
end
auxvetor=s;

auxy=y;
j=1;
for i=auxvetor
    if i~=auxvetor(length(auxvetor))
        vej=y(auxvetor(j+1),:)-centelem(ielem,:);
        % Estes condições evitam que o acos seja numero complexo.
        if (dot(vej,ve2)/(norm(vej)*norm(ve2)))>1 && abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<1e-10
            thetalef2=acos(1);
        else
            thetalef2=acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
        end
        
        vei=y(i,:)-centelem(ielem,:);
        % analiza que o K.n pertece ao primeiro quadrante
        auxquadrant1= cross(vei,ve2);
        auxquadrant2= cross(ve2,vej);
        % evita que aparição de numeros complexos
        if (dot(vei,ve2)/(norm(vei)*norm(ve2)))>1 && abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<1e-10
            % calculo do theta1
            thetalef1=acos(1);
            
        else
            thetalef1=acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
        end
        if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-10 || abs(auxquadrant2(1,3))>1e-10))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
            ksii=dot(cross(ve2,vej),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
            ksij=dot(cross(vei,ve2),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
            aux11=i;
            aux12=auxvetor(j+1);
        end
    else
        vej=y(auxvetor(1),:)-centelem(ielem,:);
        % Estes condições evitam que o acos seja numero complexo.
        if (dot(vej,ve2)/(norm(vej)*norm(ve2)))>1 && abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<1e-10
            thetalef2=acos(1);
        else
            thetalef2=acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
        end
        
        vei=y(i,:)-centelem(ielem,:);
        % analiza que o K.n pertece ao primeiro quadrante
        auxquadrant1= cross(vei,ve2);
        auxquadrant2= cross(ve2,vej);
        % evita que aparição de numeros complexos
        if (dot(vei,ve2)/(norm(vei)*norm(ve2)))>1 && abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<1e-10
            
            % calculo do theta1
            thetalef1=acos(1);
            
        else
            thetalef1=acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
        end
        if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-10 || abs(auxquadrant2(1,3))>1e-10))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
            ksii=dot(cross(ve2,vej),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
            ksij=dot(cross(vei,ve2),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
            aux11=i;
            aux12=auxvetor(1);
            
        end
        
    end
    j=j+1;
end
end
