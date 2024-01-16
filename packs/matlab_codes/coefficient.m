function [coefficient]=coefficient(F,y,kmap,numface)
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
    
    vetor1= y(ifacont,:)-centelem(lef,:);
    vetor2= Klef'*(R*(IJ/normIJ)');
    
    % Estes condições evitam que o acos seja numero complexo.
    aa0= abs(1-dot(vetor1,vetor2)/(norm(vetor1)*norm(vetor2)))<1e-10;
    bb0= abs(1-dot(vetor1,vetor2)/(norm(vetor1)*norm(vetor2)))>1e-10 ||abs(1-dot(vetor1,vetor2)/(norm(vetor1)*norm(vetor2)))==1e-10;
    
    thetalef2=aa0*acos(1)+bb0*acos(dot(vetor1,vetor2)/(norm(vetor1)*norm(vetor2)));
    %% percorrendo todos as faces dos elemento "lef"
    auxvetor=F(lef,1:klef);
    auxthetalef1=1e5;
    auxface=0;
    for i=auxvetor
        if i~=ifacont
            vetor3= y(i,:)-centelem(lef,:);
            %% analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vetor1,vetor2);
            auxquadrant2= cross(vetor2,vetor3);
            %% evita que aparição de numeros complexos
            nn0= abs(1-dot(vetor3,vetor2)/(norm(vetor3)*norm(vetor2)))<1e-10;
            mm0= abs(1-dot(vetor3,vetor2)/(norm(vetor3)*norm(vetor2)))>1e-10 ||abs(1-dot(vetor3,vetor2)/(norm(vetor3)*norm(vetor2)))==1e-10;
            %% calculo do theta1
            thetalef1=nn0*acos(1)+mm0*acos(dot(vetor3,vetor2)/(norm(vetor3)*norm(vetor2)));
             if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                     (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && (norm(auxquadrant1-auxquadrant2)>1e-10)&& ((thetalef2 + thetalef1)<pi)
%            if auxthetalef1>thetalef1 || abs(auxthetalef1-thetalef1)<1e-5    
                auxthetalef1= thetalef1;
                auxface=i;
            else
                auxthetalef1=  auxthetalef1;
                auxface=auxface;
            end
        end
    end
   
    % correspondente a face diferente do atual (em sentido anti-horario)
    coefficient(1,1,ifacont)=(norm(vetor2)*sin(thetalef2))/(norm(y(auxface,:)-centelem(lef,:))*sin(auxthetalef1+thetalef2));
    % correpondente a face atual
    coefficient(1,2,ifacont)=(norm(vetor2)*sin(auxthetalef1))/(norm(vetor1)*sin(auxthetalef1+thetalef2));
    coefficient(1,3,ifacont)=auxface;
    coefficient(1,4,ifacont)=numface(ifacont,1);
    %% satisfazendo a inedtidade
    coefficient(1,5,ifacont)= norm(((y(auxface,:)-centelem(lef,:))*coefficient(1,1,ifacont) + vetor1*coefficient(1,2,ifacont))/norm(vetor2));

end
%edgecontrib(:,:,size(bedge,1)+1:size(bedge,1)+size(inedge,1))= 0;
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
    ve1=y(iface+size(bedge,1),:)-centelem(lef,:);
    ve2=Klef*(R*(IJ/normIJ)');
    
    % Estes condições evitam que o acos seja numero complexo.
    aa=abs(1-dot(ve1,ve2)/(norm(ve1)*norm(ve2)))<1e-10;
    bb=abs(1-dot(ve1,ve2)/(norm(ve1)*norm(ve2)))>1e-10 || abs(1-dot(ve1,ve2)/(norm(ve1)*norm(ve2)))==1e-10;
    
    thetalef2=aa*acos(1)+ bb*acos( dot(ve1,ve2)/(norm(ve1)*norm(ve2)));
    %% Elemento a esquerda
    auxvetor=F(lef,1:klef);
    auxthetalef1=1e5;
    aux11=0;
    for i=auxvetor
        if i~=iface+size(bedge,1)
            
            vetor9=y(i,:)-centelem(lef,:);
            %% analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(ve1,ve2);
            auxquadrant2= cross(ve2,vetor9);
            %% evita que aparição de numeros complexos
            cc= abs(1-dot(vetor9,ve2)/(norm(vetor9)*norm(ve2)))<1e-10;
            dd= abs(1-dot(vetor9,ve2)/(norm(vetor9)*norm(ve2)))>1e-10||abs(1-dot(vetor9,ve2)/(norm(vetor9)*norm(ve2)))==1e-10;
            %% calculo do theta1
            thetalef1=cc*acos(1) + dd*acos(dot(vetor9,ve2)/(norm(vetor9)*norm(ve2)));
              if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                      (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && (norm(auxquadrant1-auxquadrant2)>1e-10)&& ((thetalef2 + thetalef1)<pi)
%            if auxthetalef1>thetalef1 || abs(auxthetalef1-thetalef1)<1e-5
                auxthetalef1=thetalef1;
                aux11=i;
            else
                auxthetalef1= auxthetalef1;
                aux11=aux11;
            end
        end
    end
    % correspondente a face anterior do atual (em sentido anti-horario)
    coefficient(1,1,iface+size(bedge,1))=(norm(ve2))*sin(thetalef2)/(norm(y(aux11,:)-centelem(lef,:))*sin(auxthetalef1+thetalef2));
    % correpondente a face atual
    coefficient(1,2,iface+size(bedge,1))=(norm(ve2))*sin(auxthetalef1)/(norm(y(iface+size(bedge,1),:)-centelem(lef,:))*sin(auxthetalef1+thetalef2));
    coefficient(1,3,iface+size(bedge,1))=aux11;
    coefficient(1,4,iface+size(bedge,1))=numface(iface+size(bedge,1),1);
    
    %% satisfazendo a inedtidade
    coefficient(1,5,iface+size(bedge,1))= norm(((y(aux11,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1)) + (y(iface+size(bedge,1),:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1)))/norm(ve2));
      
    %% Elemento a direita
    % iface+size(bedge,1): Ponto armonico na face atual
    vetor11=y(iface+size(bedge,1),:)-centelem(rel,:);
    vetor12=Krel*(R*(-IJ/normIJ)');
    
    % Estes condições evitam que o acos seja numero complexo.
    aa2=abs(1-dot(vetor11,vetor12)/(norm(vetor11)*norm(vetor12)))<1e-10;
    bb2=abs(1-dot(vetor11,vetor12)/(norm(vetor11)*norm(vetor12)))>1e-10 || abs(1-dot(vetor11,vetor12)/(norm(vetor11)*norm(vetor12)))==1e-10;
    
    thetarel2=aa2*acos(1)+bb2*acos( dot(vetor11,vetor12)/(norm(vetor11)*norm(vetor12)));
    
    auxvetorel=F(rel,1:krel);
    auxthetarel1=1e5;
    aux22=0;
    for ii=auxvetorel
        if ii~=iface+size(bedge,1)
            
            vetor15=y(ii,:)-centelem(rel,:);
            %% analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vetor11,vetor12);
            auxquadrant2= cross(vetor12,vetor15);
            %% evita que aparição de numeros complexos
            cc2=abs(1-dot(vetor15,vetor12)/(norm(vetor15)*norm(vetor12)))<1e-10;
            dd2=abs(1-dot(vetor15,vetor12)/(norm(vetor15)*norm(vetor12)))>1e-10 || abs(1-dot(vetor15,vetor12)/(norm(vetor15)*norm(vetor12)))==1e-10;
            %% calculo do theta1
            thetarel1=cc2*acos(1)+dd2*acos( dot(vetor15,vetor12)/(norm(vetor15)*norm(vetor12)));
            
              if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>1e-5 || abs(auxquadrant2(1,3))>1e-5))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                      (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && (norm(auxquadrant1-auxquadrant2)>1e-10) && ((thetarel2 + thetarel1)<pi)
%            if auxthetarel1>thetarel1 || abs(auxthetarel1-thetarel1)<1e-5
                auxthetarel1=thetarel1;
                aux22=ii;
            else
                auxthetarel1=auxthetarel1;
                aux22=aux22;
            end
        end
    end
    % correspondente a face anterior do atual (em sentido anti-horario)
    coefficient(2,1,iface+size(bedge,1))=(norm(vetor12))*sin(thetarel2)/(norm(y(aux22,:)-centelem(rel,:))*sin(auxthetarel1+thetarel2));
    % correpondente a face atual
    coefficient(2,2,iface+size(bedge,1))=(norm(vetor12))*sin(auxthetarel1)/(norm(y(iface+size(bedge,1),:)-centelem(rel,:))*sin(auxthetarel1+thetarel2));
    coefficient(2,3,iface+size(bedge,1))=aux22;
    coefficient(2,4,iface+size(bedge,1))=numface(iface+size(bedge,1),1);
    %% satisfazendo a inedtidade
    coefficient(2,5,iface+size(bedge,1))=norm(((y(iface+size(bedge,1),:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))+ (y(aux22,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1)))/norm(vetor12));
  
    
end

end