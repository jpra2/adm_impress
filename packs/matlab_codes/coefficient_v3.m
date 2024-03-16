function [coefficient]=coefficient_v3(F,y,kmap,numface)
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
    auxvetor=F(lef,1:klef);
    j=1;
    for i=auxvetor
        if i~=auxvetor(length(auxvetor))
            vej=y(auxvetor(j+1),:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            aa=abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<1e-10;
            bb=abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))>1e-10 || abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))==1e-10;
            
            thetalef2=aa*acos(1)+ bb*acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
            
            vei=y(i,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vej,ve2);
            auxquadrant2= cross(ve2,vei);
            % evita que aparição de numeros complexos
            cc= abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<1e-10;
            dd= abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))>1e-10||abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))==1e-10;
            % calculo do theta1
            thetalef1=cc*acos(1) + dd*acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                auxthetalef1=thetalef1;
                auxthetalef2=thetalef2;
                aux11=i;
                aux12=auxvetor(j+1);

            end
        else
            vej=y(auxvetor(1),:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            aa=abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<1e-10;
            bb=abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))>1e-10 || abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))==1e-10;
            
            thetalef2=aa*acos(1)+ bb*acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
            
            vei=y(i,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vej,ve2);
            auxquadrant2= cross(ve2,vei);
            % evita que aparição de numeros complexos
            cc= abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<1e-10;
            dd= abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))>1e-10||abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))==1e-10;
            % calculo do theta1
            thetalef1=cc*acos(1) + dd*acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                auxthetalef1=thetalef1;
                auxthetalef2=thetalef2;
                aux11=i;
                aux12=auxvetor(1);
               
            end
            
        end
        j=j+1;
    end
    a=[aux11 aux12];
    c=find(a==ifacont);
    if c~=0
        if a(c)==aux11
            % atribuindo valores a os coeficientes
            coefficient(1,1,ifacont)=(norm(ve2))*sin(auxthetalef2)/(norm(y(aux11,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
            coefficient(1,2,ifacont)=(norm(ve2))*sin(auxthetalef1)/(norm(y(aux12,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
            coefficient(1,3,ifacont)=0;
            % indexando as faces respetivamente
            coefficient(1,4,ifacont)=aux11;
            coefficient(1,5,ifacont)=aux12;
            coefficient(1,6,ifacont)=0;
            % verificando a identidade
            coefficient(1,7,ifacont)=norm(((y(aux11,:)-centelem(lef,:))*coefficient(1,1,ifacont)+ (y(aux12,:)-centelem(lef,:))*coefficient(1,2,ifacont))/norm(ve2));
            if abs(coefficient(1,7,ifacont)-1)<1e-10
            else
                disp('error')
            end
        else
            % atribuindo valores a os coeficientes
            coefficient(1,1,ifacont)=(norm(ve2))*sin(auxthetalef1)/(norm(y(aux12,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
            coefficient(1,2,ifacont)=(norm(ve2))*sin(auxthetalef2)/(norm(y(aux11,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
            coefficient(1,3,ifacont)=0;
            % indexando as faces respetivamente
            coefficient(1,4,ifacont)=aux12;
            coefficient(1,5,ifacont)=aux11;
            coefficient(1,6,ifacont)=0;
            % verificando a identidade
            coefficient(1,7,ifacont)=norm(((y(aux12,:)-centelem(lef,:))*coefficient(1,1,ifacont)+ (y(aux11,:)-centelem(lef,:))*coefficient(1,2,ifacont))/norm(ve2));
            if abs(coefficient(1,7,ifacont)-1)<1e-10
            else
                disp('error')
            end
        end
        
    else
        % atribuindo valores a os coeficientes
        coefficient(1,1,ifacont)=0;
        % correpondente a face atual
        coefficient(1,2,ifacont)=(norm(ve2))*sin(auxthetalef1)/(norm(y(aux12,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
        coefficient(1,3,ifacont)=(norm(ve2))*sin(auxthetalef2)/(norm(y(aux11,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
        % indexando as faces respetivamente
        coefficient(1,4,ifacont)=0;
        coefficient(1,5,ifacont)=aux12;
        coefficient(1,6,ifacont)=aux11;
        % verificando a identidade
        coefficient(1,7,ifacont)=norm(((y(aux12,:)-centelem(lef,:))*coefficient(1,2,ifacont)+ (y(aux11,:)-centelem(lef,:))*coefficient(1,3,ifacont))/norm(ve2));
        if abs(coefficient(1,7,ifacont)-1)<1e-10
        else
            disp('error')
        end
    end
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
    % ponto armonico no face interior atual
    % primeira maneira de calcular kn
    %     auxvetor=F(lef,1:klef);
    %     c=find(auxvetor==iface+size(bedge,1));
    %     if c==1
    %         d=auxvetor(length(auxvetor));
    %     else
    %        d=auxvetor(c-1);
    %     end
    %     [ve2]=calconormal(inedge(iface,1),inedge(iface,2),auxvetor(c), d,y,Klef,lef);
    % sugunda maneira de calcular o kn
    %       Knlef=((R*IJ')'*Klef*(R*IJ'))/norm(IJ)^2;
    %       Ktlef=((R*IJ')'*Klef*IJ')/norm(IJ)^2;
    %
    %     ve2=Ktlef*(IJ'/norm(IJ))+Knlef*(R*(IJ/normIJ)');
    
    % terceira maneira de calcular o kn
    ve2=Klef*(R*(IJ/normIJ)');
    % faces na que compoem o elemento    
    auxvetor=F(lef,1:klef);
    
    j=1;
    for i=auxvetor
        if i~=auxvetor(length(auxvetor))
            vej=y(auxvetor(j+1),:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            aa=abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<1e-10;
            bb=abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))>1e-10 || abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))==1e-10;
            
            thetalef2=aa*acos(1)+ bb*acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
            
            vei=y(i,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vej,ve2);
            auxquadrant2= cross(ve2,vei);
            % evita que aparição de numeros complexos
            cc= abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<1e-10;
            dd= abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))>1e-10||abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))==1e-10;
            % calculo do theta1
            thetalef1=cc*acos(1) + dd*acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                auxthetalef1=thetalef1;
                auxthetalef2=thetalef2;
                aux11=i;
                aux12=auxvetor(j+1);
                
                
            end
        else
            vej=y(auxvetor(1),:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            aa=abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<1e-10;
            bb=abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))>1e-10 || abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))==1e-10;
            
            thetalef2=aa*acos(1)+ bb*acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
            
            vei=y(i,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vej,ve2);
            auxquadrant2= cross(ve2,vei);
            % evita que aparição de numeros complexos
            cc= abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<1e-10;
            dd= abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))>1e-10||abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))==1e-10;
            % calculo do theta1
            thetalef1=cc*acos(1) + dd*acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                auxthetalef1=thetalef1;
                auxthetalef2=thetalef2;
                aux11=i;
                aux12=auxvetor(1);
               
            end
            
        end
        j=j+1;
    end
    a=[aux11 aux12];
    c=find(a==iface+size(bedge,1));
    if c~=0
        if a(c)==aux11
            % atribuindo valores a os coeficientes
            coefficient(1,1,iface+size(bedge,1))=(norm(ve2))*sin(auxthetalef2)/(norm(y(aux11,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
            coefficient(1,2,iface+size(bedge,1))=(norm(ve2))*sin(auxthetalef1)/(norm(y(aux12,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
            coefficient(1,3,iface+size(bedge,1))=0;
            % indexando as faces respetivamente
            coefficient(1,4,iface+size(bedge,1))=aux11;
            coefficient(1,5,iface+size(bedge,1))=aux12;
            coefficient(1,6,iface+size(bedge,1))=0;
            % verificando a identidade
            coefficient(1,7,iface+size(bedge,1))=norm(((y(aux11,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1))+ (y(aux12,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1)))/norm(ve2));
            if abs(coefficient(1,7,iface+size(bedge,1))-1)<1e-10
            else
                disp('error')
            end
        else
            % atribuindo valores a os coeficientes
            coefficient(1,1,iface+size(bedge,1))=(norm(ve2))*sin(auxthetalef1)/(norm(y(aux12,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
            coefficient(1,2,iface+size(bedge,1))=(norm(ve2))*sin(auxthetalef2)/(norm(y(aux11,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
            coefficient(1,3,iface+size(bedge,1))=0;
            % indexando as faces respetivamente
            coefficient(1,4,iface+size(bedge,1))=aux12;
            coefficient(1,5,iface+size(bedge,1))=aux11;
            coefficient(1,6,iface+size(bedge,1))=0;
            % verificando a identidade
            coefficient(1,7,iface+size(bedge,1))=norm(((y(aux12,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1))+ (y(aux11,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1)))/norm(ve2));
            if abs(coefficient(1,7,iface+size(bedge,1))-1)<1e-10
            else
                disp('error')
            end
        end
        
    else
        % atribuindo valores a os coeficientes
        coefficient(1,1,iface+size(bedge,1))=0;
        % correpondente a face atual
        coefficient(1,2,iface+size(bedge,1))=(norm(ve2))*sin(auxthetalef1)/(norm(y(aux12,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
        coefficient(1,3,iface+size(bedge,1))=(norm(ve2))*sin(auxthetalef2)/(norm(y(aux11,:)-centelem(lef,:))*sin(auxthetalef1+auxthetalef2));
        % indexando as faces respetivamente
        coefficient(1,4,iface+size(bedge,1))=0;
        coefficient(1,5,iface+size(bedge,1))=aux12;
        coefficient(1,6,iface+size(bedge,1))=aux11;
        % verificando a identidade
        coefficient(1,7,iface+size(bedge,1))=norm(((y(aux12,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1))+ (y(aux11,:)-centelem(lef,:))*coefficient(1,3,iface+size(bedge,1)))/norm(ve2));
        if abs(coefficient(1,7,iface+size(bedge,1))-1)<1e-10
        else
            disp('error')
        end
    end
    clear aux12 aux11
    
    %% Elemento a direita
    % iface+size(bedge,1): Ponto armonico na face atual
    %    vetor11=y(iface+size(bedge,1),:)-centelem(rel,:);
    % primeira de calcular os kn
    %     c=find(auxvetor==iface+size(bedge,1));
    %     if c==1
    %         d=auxvetor(length(auxvetor));
    %     else
    %        d=auxvetor(c-1);
    %     end
    %     [vetor12]=calconormal(inedge(iface,2),inedge(iface,1),auxvetor(c), d,y,Krel,rel);
    % segunda maneira de calcular os kn
    % Knrel=((R*IJ')'*Krel*(R*IJ'))/norm(IJ)^2;
    % Ktrel=((R*IJ')'*Krel*IJ')/norm(IJ)^2;
    %
    % vetor12=Ktrel*(-IJ'/norm(IJ))+Knrel*(R*(-IJ/normIJ)');
    % terceira maneira de calcular os kn
    vetor12=Krel*(R*(-IJ/normIJ)');
    auxvetorel=F(rel,1:krel);
    
    j=1;
    for ii=auxvetorel
        if ii~=auxvetorel(length(auxvetorel))
            vetorj=y(auxvetorel(j+1),:)-centelem(rel,:);
            % Estes condições evitam que o acos seja numero complexo.
            aa2=abs(1-dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))<1e-10;
            bb2=abs(1-dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))>1e-10 || abs(1-dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))==1e-10;
            % calculo do theta2
            thetarel2=aa2*acos(1)+bb2*acos( dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)));

            vetori=y(ii,:)-centelem(rel,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vetorj,vetor12);
            auxquadrant2= cross(vetor12,vetori);
            % evita que aparição de numeros complexos
            cc2=abs(1-dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))<1e-10;
            dd2=abs(1-dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))>1e-10 || abs(1-dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))==1e-10;
            % calculo do theta1
            thetarel1=cc2*acos(1)+dd2*acos( dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)));
            
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0))  && ((thetarel2 + thetarel1)<pi)
                auxthetarel1=thetarel1;
                auxthetarel2=thetarel2;
                auxi=ii;
                auxj=auxvetorel(j+1);
                
                
            end
            
        else
            vetorj=y(auxvetorel(1),:)-centelem(rel,:);
            % Estes condições evitam que o acos seja numero complexo.
            aa2=abs(1-dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))<1e-10;
            bb2=abs(1-dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))>1e-10 || abs(1-dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))==1e-10;
            % calculo do theta2
            thetarel2=aa2*acos(1)+bb2*acos( dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)));
            
            vetori=y(ii,:)-centelem(rel,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vetorj,vetor12);
            auxquadrant2= cross(vetor12,vetori);
            % evita que aparição de numeros complexos
            cc2=abs(1-dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))<1e-10;
            dd2=abs(1-dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))>1e-10 || abs(1-dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))==1e-10;
            % calculo do theta1
            thetarel1=cc2*acos(1)+dd2*acos( dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)));
            
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0))  && ((thetarel2 + thetarel1)<pi)
                auxthetarel1=thetarel1;
                auxthetarel2=thetarel2;
                auxi=ii;
                auxj=auxvetorel(1);
                
                
            end
        end
        j=j+1;
        
    end
    a=[auxi auxj];
    c=find(a==iface+size(bedge,1));
    if c~=0
        if a(c)==auxi
            % atribuindo valores a os coeficientes
            coefficient(2,1,iface+size(bedge,1))=(norm(vetor12))*sin(auxthetarel2)/(norm(y(auxi,:)-centelem(rel,:))*sin(auxthetarel1+auxthetarel2));
            coefficient(2,2,iface+size(bedge,1))=(norm(vetor12))*sin(auxthetarel1)/(norm(y(auxj,:)-centelem(rel,:))*sin(auxthetarel1+auxthetarel2));
            coefficient(2,3,iface+size(bedge,1))=0;
            % indexando as faces respetivamente
            coefficient(2,4,iface+size(bedge,1))=auxi;
            coefficient(2,5,iface+size(bedge,1))=auxj;
            coefficient(2,6,iface+size(bedge,1))=0;
            % verificando a identidade
            coefficient(2,7,iface+size(bedge,1))=norm(((y(auxi,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1))+ (y(auxj,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1)))/norm(vetor12));
            if abs(coefficient(2,7,iface+size(bedge,1))-1)<1e-10
            else
               
                disp('error')
                
            end
            
        else
            % atribuindo valores a os coeficientes
            coefficient(2,1,iface+size(bedge,1))=(norm(vetor12))*sin(auxthetarel1)/(norm(y(auxj,:)-centelem(rel,:))*sin(auxthetarel1+auxthetarel2));
            coefficient(2,2,iface+size(bedge,1))=(norm(vetor12))*sin(auxthetarel2)/(norm(y(auxi,:)-centelem(rel,:))*sin(auxthetarel1+auxthetarel2));
            coefficient(2,3,iface+size(bedge,1))=0;
            % indexando as faces respetivamente
            coefficient(2,4,iface+size(bedge,1))=auxj;
            coefficient(2,5,iface+size(bedge,1))=auxi;
            coefficient(2,6,iface+size(bedge,1))=0;
            % verificando a identidade
            coefficient(2,7,iface+size(bedge,1))=norm(((y(auxi,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))+ (y(auxj,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1)))/norm(vetor12));
            if abs(coefficient(2,7,iface+size(bedge,1))-1)<1e-10
            else
                disp('error')
            end 
        end 
    else
        % atribuindo valores a os coeficientes
        coefficient(2,1,iface+size(bedge,1))=0;
        coefficient(2,2,iface+size(bedge,1))=(norm(vetor12))*sin(auxthetarel1)/(norm(y(auxj,:)-centelem(rel,:))*sin(auxthetarel1+auxthetarel2));
        coefficient(2,3,iface+size(bedge,1))=(norm(vetor12))*sin(auxthetarel2)/(norm(y(auxi,:)-centelem(rel,:))*sin(auxthetarel1+auxthetarel2));
        % indexando as faces respetivamente
        coefficient(2,4,iface+size(bedge,1))=0;
        coefficient(2,5,iface+size(bedge,1))=auxj;
        coefficient(2,6,iface+size(bedge,1))=auxi;
        % verificando a inedtidade
        coefficient(2,7,iface+size(bedge,1))=norm(((y(auxj,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))+ (y(auxi,:)-centelem(rel,:))*coefficient(2,3,iface+size(bedge,1)))/norm(vetor12));
        if abs(coefficient(2,7,iface+size(bedge,1))-1)<1e-10
        else
            disp('error')
        end
    end
    clear auxj auxi
       
end

end