function [coeff]=coefficient_v1(F,y,kmap)
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
 
     
    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);

    %% ej(iface) do elemento a esquerda
    a= Klef'*(R*(IJ/normIJ)');
    
    b= y(ifacont,:)-centelem(lef,:);
    %% percorrendo todos as faces dos elemento "lef"
    auxvetor=F(lef,1:klef);
    for i=auxvetor
        if i~=ifacont
            c=y(i,:)-centelem(lef,:);
            aux1=cross(c,b);
            aux2=cross(a',b);
            aux3=cross(c,a');
            
            Csii=dot(aux2,[0,0,sign(aux1(1,3))*1])/dot(aux1,[0,0,sign(aux1(1,3))*1]);
            
            Csij=dot(aux3,[0,0,sign(aux1(1,3))*1])/dot(aux1,[0,0,sign(aux1(1,3))*1]);
            
            r= Csii*c + Csij*b;
            
            if (norm(r-a')<1e-10) && (Csii>0 ||Csii==0) && (Csij>0 ||Csij==0)
                
                % correspondente a face diferente do atual
                coeff(1,1,ifacont)=Csii;
                % correpondente a face atual
                coeff(1,2,ifacont)=Csij;
                coeff(1,3,ifacont)=i;
                coeff(1,4,ifacont)=norm(r-a');
                break
                
            end
        end
    end
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
    
   

% %Essa é UMA maneira de construir os tensores.

    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);

    % tensor a elemento a direita

    Krel(1,1)=kmap(elem(rel,5),2);
    Krel(1,2)=kmap(elem(rel,5),3);
    Krel(2,1)=kmap(elem(rel,5),4);
    Krel(2,2)=kmap(elem(rel,5),5);
    %% ej(iface) do elemento a esquerda
    % ponto armonico no face interior atual
    alef=Klef'*(R*(IJ/normIJ)');
    
    blef=y(iface+size(bedge,1),:)-centelem(lef,:);
    
    %% Elemento a esquerda
    auxvetorlef=F(lef,1:klef);
    
    for i=auxvetorlef
        if i~=iface+size(bedge,1)
            
            clef=y(i,:)-centelem(lef,:);
            aux=cross(clef,blef);
            Csiilef=dot(cross(alef',blef),[0,0,sign(aux(1,3))*1])/dot(aux,[0,0,sign(aux(1,3))*1]);
            
            Csijlef=dot(cross(clef,alef'),[0,0,sign(aux(1,3))*1])/dot(aux,[0,0,sign(aux(1,3))*1]);
            
            rlef=Csiilef*clef + Csijlef*blef;
            
            if (norm(alef'-rlef)<1e-10) && (Csiilef>0 || Csiilef ==0) && (Csijlef>0 || Csijlef==0)
                % correspondente a face anterior do atual
                coeff(1,1,iface+size(bedge,1))=Csiilef;
                % correpondente a face atual
                coeff(1,2,iface+size(bedge,1))=Csijlef;
                coeff(1,3,iface+size(bedge,1))=i;
                coeff(1,4,iface+size(bedge,1))=norm(rlef-alef'); % verifica quando o vetor é unitario
                break
            end
        end
    end
    
    %% Elemento a direita
    % iface+size(bedge,1): Ponto armonico na face atual
    
    arel=Krel'*(R*(-IJ/normIJ)');
    
    brel=y(iface+size(bedge,1),:)-centelem(rel,:);
    
    auxvetorel=F(rel,1:krel);
    
    for ii=auxvetorel
        if ii~=iface+size(bedge,1)
            
            crel=y(ii,:)-centelem(rel,:);
            aux=cross(crel,brel);
            Csiirel=dot(cross(arel,brel),[0,0,sign(aux(1,3))*1])/dot(aux,[0,0,sign(aux(1,3))*1]);
            
            Csijrel=dot(cross(crel,arel),[0,0,sign(aux(1,3))*1])/dot(aux,[0,0,sign(aux(1,3))*1]);
            
            rrel=Csiirel*crel + Csijrel*brel;
            if (norm(rrel-arel')<1e-10) && (Csiirel>0 ||Csiirel==0) && (Csijrel>0 || Csijrel==0)
                % correpondente a face diferente do atual
                coeff(2,1,iface+size(bedge,1))=Csiirel;
                % correspondente a face atual
                coeff(2,2,iface+size(bedge,1))=Csijrel;
                coeff(2,3,iface+size(bedge,1))=ii;
                coeff(2,4,iface+size(bedge,1))=norm(rrel-arel'); % verifica quando o vetor é unitario
                break
            end
        end
    end
end

end