function [coefficient,calnormface]=coefficientLPSangle(kmap)
global inedge bedge coord elem centelem
Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0;0 0 0];

for ifacont=1:size(bedge,1)
    
    lef=bedge(ifacont,3);
    klef=elementype(lef);
    
    IJ=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    normIJ=norm(IJ);
    %Essa é UMA maneira de construir os tensores.
    
    Klef(1,1)=kmap(elem(bedge(ifacont,3),5),2);
    Klef(1,2)=kmap(elem(bedge(ifacont,3),5),3);
    Klef(2,1)=kmap(elem(bedge(ifacont,3),5),4);
    Klef(2,2)=kmap(elem(bedge(ifacont,3),5),5);
    
    %% ej(iface) do elemento a esquerda
    
    ve2= Klef'*(R*(IJ/normIJ)');
    %% percorrendo todos as faces dos elemento "lef"
    auxvetor=elem(lef,1:klef);
    j=1;
    for i=auxvetor
        if j==length(auxvetor)
            no1=i;
            no2=auxvetor(1);
            vj=coord(no2,:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,ve2)/(norm(vj)*norm(ve2))>1 && abs(1-dot(vj,ve2)/(norm(vj)*norm(ve2)))<1e-10
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vj,ve2)/(norm(vj)*norm(ve2)));
            end
            vi=coord(no1,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vj,ve2);
            auxquadrant2= cross(ve2,vi);
            % evita que aparição de numeros complexos
            if dot(vi,ve2)/(norm(vi)*norm(ve2))>1 && abs(1-dot(vi,ve2)/(norm(vi)*norm(ve2)))<1e-10
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vi,ve2)/(norm(vi)*norm(ve2)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-16 || abs(auxquadrant2(1,3))>1e-16))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                
                aux11=no1;
                aux12=no2;
                % atribuindo valores a os coeficientes
                coefficient(1,1,ifacont)=(norm(ve2))*sin(thetalef2)/(norm(coord(aux11,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                coefficient(1,2,ifacont)=(norm(ve2))*sin(thetalef1)/(norm(coord(aux12,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                
                % indexando as faces respetivamente
                coefficient(1,3,ifacont)=aux11;
                coefficient(1,4,ifacont)=aux12;
                
                % verificando a identidade
                coefficient(1,5,ifacont)=norm((coord(aux11,:)-centelem(lef,:))*coefficient(1,1,ifacont)+...
                    (coord(aux12,:)-centelem(lef,:))*coefficient(1,2,ifacont)-ve2');
                if abs(coefficient(1,5,ifacont))<1e-8
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
                
            end
        else
            no1=i;
            no2=auxvetor(j+1);
            vj= coord(no2,:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,ve2)/(norm(vj)*norm(ve2))>1 && abs(1-dot(vj,ve2)/(norm(vj)*norm(ve2)))<1e-10
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vj,ve2)/(norm(vj)*norm(ve2)));
            end
            vi= coord(no1,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vi,ve2);
            auxquadrant2= cross(ve2,vj);
            % evita que aparição de numeros complexos
            if dot(vi,ve2)/(norm(vi)*norm(ve2))>1 && abs(1-dot(vi,ve2)/(norm(vi)*norm(ve2)))<1e-10
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vi,ve2)/(norm(vi)*norm(ve2)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-16 || abs(auxquadrant2(1,3))>1e-16))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                
                aux11=no1;
                aux12=no2;
                % atribuindo valores a os coeficientes
                coefficient(1,1,ifacont)=(norm(ve2))*sin(thetalef2)/(norm(coord(aux11,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                coefficient(1,2,ifacont)=(norm(ve2))*sin(thetalef1)/(norm(coord(aux12,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                
                % indexando as faces respetivamente
                coefficient(1,3,ifacont)=aux11;
                coefficient(1,4,ifacont)=aux12;
                
                % verificando a identidade
                coefficient(1,5,ifacont)=norm((coord(aux11,:)-centelem(lef,:))*coefficient(1,1,ifacont)+(coord(aux12,:)-centelem(lef,:))*coefficient(1,2,ifacont)-ve2');
                if abs(coefficient(1,5,ifacont))<1e-8
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
            end
            
        end
        j=j+1;
    end
 calnormface(ifacont,1)=normIJ; 
    clear aux12 aux11
end

%% Faces interiores
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    klef=elementype(lef);
    krel=elementype(rel);
    
    IJ=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    normIJ=norm(coord(inedge(iface,2),:)-coord(inedge(iface,1),:));
    %Essa é UMA maneira de construir os tensores.
    calnormface(iface,1)=normIJ; 
    Klef(1,1)=kmap(elem(inedge(iface,3),5),2);
    Klef(1,2)=kmap(elem(inedge(iface,3),5),3);
    Klef(2,1)=kmap(elem(inedge(iface,3),5),4);
    Klef(2,2)=kmap(elem(inedge(iface,3),5),5);
    
    %tensor a elemento a direita
    
    Krel(1,1)=kmap(elem(inedge(iface,4),5),2);
    Krel(1,2)=kmap(elem(inedge(iface,4),5),3);
    Krel(2,1)=kmap(elem(inedge(iface,4),5),4);
    Krel(2,2)=kmap(elem(inedge(iface,4),5),5);
    
    %% ej(iface) do elemento a esquerda
    
    % terceira maneira de calcular o kn
    ve2=Klef*(R*(IJ/normIJ)');
    % faces na que compoem o elemento
    auxvetor=elem(lef,1:klef);
    
    j=1;
    for i=auxvetor
        if j==length(auxvetor)
            no1=i;
            no2=auxvetor(1);
            vj=coord(no2,:)-centelem(lef,:);
            
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,ve2)/(norm(vj)*norm(ve2))>1 && abs(1-dot(vj,ve2)/(norm(vj)*norm(ve2)))<1e-10
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vj,ve2)/(norm(vj)*norm(ve2)));
            end
            vi=coord(no1,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vj,ve2);
            auxquadrant2= cross(ve2,vi);
            % evita que aparição de numeros complexos
            if dot(vi,ve2)/(norm(vi)*norm(ve2))>1 && abs(1-dot(vi,ve2)/(norm(vi)*norm(ve2)))<1e-10
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vi,ve2)/(norm(vi)*norm(ve2)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-16 || abs(auxquadrant2(1,3))>1e-16))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                
                aux11=no1;
                aux12=no2;
                % atribuindo valores a os coeficientes
                coefficient(1,1,iface+size(bedge,1))=(norm(ve2))*sin(thetalef2)/(norm(coord(aux11,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                coefficient(1,2,iface+size(bedge,1))=(norm(ve2))*sin(thetalef1)/(norm(coord(aux12,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                
                % indexando as faces respetivamente
                coefficient(1,3,iface+size(bedge,1))=aux11;
                coefficient(1,4,iface+size(bedge,1))=aux12;
                
                % verificando a identidade
                coefficient(1,5,iface+size(bedge,1))=norm((coord(aux11,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1))+...
                    (coord(aux12,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1))- ve2');
                if abs(coefficient(1,5,iface+size(bedge,1)))<1e-8
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
                
            end
        else
            no1=i;
            no2=auxvetor(j+1);
            vj= coord(no2,:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,ve2)/(norm(vj)*norm(ve2))>1 && abs(1-dot(vj,ve2)/(norm(vj)*norm(ve2)))<1e-10
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vj,ve2)/(norm(vj)*norm(ve2)));
            end
            vi= coord(no1,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vj,ve2);
            auxquadrant2= cross(ve2,vi);
            % evita que aparição de numeros complexos
            if dot(vi,ve2)/(norm(vi)*norm(ve2))>1 && abs(1-dot(vi,ve2)/(norm(vi)*norm(ve2)))<1e-10
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vi,ve2)/(norm(vi)*norm(ve2)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-16 || abs(auxquadrant2(1,3))>1e-16))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                
                aux11=no1;
                aux12=no2;
                % atribuindo valores a os coeficientes
                coefficient(1,1,iface+size(bedge,1))=(norm(ve2))*sin(thetalef2)/(norm(coord(aux11,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                coefficient(1,2,iface+size(bedge,1))=(norm(ve2))*sin(thetalef1)/(norm(coord(aux12,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                
                % indexando as faces respetivamente
                coefficient(1,3,iface+size(bedge,1))=aux11;
                coefficient(1,4,iface+size(bedge,1))=aux12;
                
                % verificando a identidade
                coefficient(1,5,iface+size(bedge,1))=norm((coord(aux11,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1))+...
                    (coord(aux12,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1))- ve2');
                if abs(coefficient(1,5,iface+size(bedge,1)))<1e-8
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
            end
            
        end
        j=j+1;
    end
    
    clear aux12 aux11
    
    %% Elemento a direita
    
    vetor12=Krel*(R*(-IJ/normIJ)');
    auxvetor=elem(rel,1:krel);
    
    j=1;
    for ii=auxvetor
        if j==length(auxvetor)
            no1=ii;
            no2=auxvetor(1);
            vj=coord(no2,:)-centelem(rel,:);
            
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,vetor12)/(norm(vj)*norm(vetor12))>1 && abs(1-dot(vj,vetor12)/(norm(vj)*norm(vetor12)))<1e-10
                % calculo do theta2
                thetarel2=acos(1);
            else
                thetarel2=acos( dot(vj,vetor12)/(norm(vj)*norm(vetor12)));
            end
            vi=coord(no1,:)-centelem(rel,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vi,vetor12);
            auxquadrant2= cross(vetor12,vj);
            % evita que aparição de numeros complexos
            if dot(vi,vetor12)/(norm(vi)*norm(vetor12))>1 && abs(1-dot(vi,vetor12)/(norm(vi)*norm(vetor12)))<1e-10
                % calculo do theta1
                thetarel1=acos(1);
            else
                thetarel1=acos( dot(vi,vetor12)/(norm(vi)*norm(vetor12)));
                
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-16 || abs(auxquadrant2(1,3))>1e-16))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0))  && ((thetarel2 + thetarel1)<pi)
                
                auxi=no1;
                auxj=no2;
                % atribuindo valores a os coeficientes
                coefficient(2,1,iface+size(bedge,1))=(norm(vetor12))*sin(thetarel2)/(norm(coord(auxi,:)-centelem(rel,:))*sin(thetarel1+thetarel2));
                coefficient(2,2,iface+size(bedge,1))=(norm(vetor12))*sin(thetarel1)/(norm(coord(auxj,:)-centelem(rel,:))*sin(thetarel1+thetarel2));
                
                % indexando as faces respetivamente
                coefficient(2,3,iface+size(bedge,1))=auxi;
                coefficient(2,4,iface+size(bedge,1))=auxj;
                
                % verificando a identidade
                coefficient(2,5,iface+size(bedge,1))=norm((coord(auxi,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1))+...
                    (coord(auxj,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))- vetor12');
                if abs(coefficient(2,5,iface+size(bedge,1)))<1e-8
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                    
                end
                
            end
            
        else
            no1=ii;
            no2=auxvetor(j+1);
            vj=coord(no2,:)-centelem(rel,:);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,vetor12)/(norm(vj)*norm(vetor12))>1 && abs(1-dot(vj,vetor12)/(norm(vj)*norm(vetor12)))<1e-10
                % calculo do theta2
                thetarel2=acos(1);
            else
                thetarel2=acos( dot(vj,vetor12)/(norm(vj)*norm(vetor12)));
            end
            
            vi= coord(no1,:)-centelem(rel,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vi,vetor12);
            auxquadrant2= cross(vetor12,vj);
            % evita que aparição de numeros complexos
            if dot(vi,vetor12)/(norm(vi)*norm(vetor12))>1 && abs(1-dot(vi,vetor12)/(norm(vi)*norm(vetor12)))<1e-10
                % calculo do theta1
                thetarel1=acos(1);
            else
                thetarel1=acos( dot(vi,vetor12)/(norm(vi)*norm(vetor12)));
                
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-16 || abs(auxquadrant2(1,3))>1e-16))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0))  && ((thetarel2 + thetarel1)<pi)
                
                auxi=no1;
                auxj=no2;
                % atribuindo valores a os coeficientes
                coefficient(2,1,iface+size(bedge,1))=(norm(vetor12))*sin(thetarel2)/(norm(coord(auxi,:)-centelem(rel,:))*sin(thetarel1+thetarel2));
                coefficient(2,2,iface+size(bedge,1))=(norm(vetor12))*sin(thetarel1)/(norm(coord(auxj,:)-centelem(rel,:))*sin(thetarel1+thetarel2));
                
                % indexando as faces respetivamente
                coefficient(2,3,iface+size(bedge,1))=auxi;
                coefficient(2,4,iface+size(bedge,1))=auxj;
                
                % verificando a identidade
                coefficient(2,5,iface+size(bedge,1))=norm((coord(auxi,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1))+...
                    (coord(auxj,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))- vetor12');
                if abs(coefficient(2,5,iface+size(bedge,1)))<1e-8
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
                
            end
        end
        j=j+1;
        
    end
    
    clear auxj auxi
    
end

end