function [G,g]=gravitation(kmap,gravelem,gravface)
global inedge bedge elem centelem coord normals strategy
Klef=zeros(3,3);
Krel=zeros(3,3);
G=zeros(size(elem,1),1);
R=[0 1 0;-1 0 0 ; 0 0 0];
K1=zeros(2,2);
K2=zeros(2,2);
K=zeros(2,2);
for ifacont=1:size(bedge,1)
    % elemento a esquerda
    lef=bedge(ifacont,3);
    % vetor de orientacao da face em questao
    ve1=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    
    % tensor de permeabilidade do elemento a esquerda
    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);
    
    Keq=Klef;
    if  strcmp(strategy,'starnoni')
        if bedge(ifacont,5)>200
            g(ifacont,1)=dot((R*ve1')'*Keq,(gravface(ifacont,:)));
        else
            g(ifacont,1)=dot((R*ve1')'*Keq,(gravelem(lef,:)));
        end
    else
        g(ifacont,1)=dot((R*ve1')'*Keq,(gravface(ifacont,:)));
    end
    G(lef,1)=G(lef,1)-g(ifacont,1);
    
end
for iface=1:size(inedge,1)
    % elementos a esquerda e a direita
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    
    % calculo do ponto meio da face
    vm=(coord(inedge(iface,1),:)+coord(inedge(iface,2),:))*0.5;
    
    % calculo da distancia do centro ao ponto meio da face
    dj1=norm(centelem(lef,:)-vm);
    dj2=norm(centelem(rel,:)-vm);
    
    
    if  strcmp(strategy,'starnoni')
        % tensor de permeabilidade do elemento a esquerda
        R1=[0 1 ;-1 0 ];
        K1(1,1)=kmap(elem(lef,5),2);
        K1(1,2)=kmap(elem(lef,5),3);
        K1(2,1)=kmap(elem(lef,5),4);
        K1(2,2)=kmap(elem(lef,5),5);
        
        % tensor de permeabilidade do elemento a direita
        
        K2(1,1)=kmap(elem(rel,5),2);
        K2(1,2)=kmap(elem(rel,5),3);
        K2(2,1)=kmap(elem(rel,5),4);
        K2(2,2)=kmap(elem(rel,5),5);
        vd1=coord(inedge(iface,2),1:2)-coord(inedge(iface,1),1:2);
        Keq=inv((dj1*inv(K1)+dj2*inv(K2))); % equation 21
        graveq=((dj1*gravelem(lef,1:2)+dj2*gravelem(rel,1:2))'); % equation 22
        g(iface+size(bedge,1),1)=dot(((R1*vd1')')*Keq, graveq);% equation 20
    else
        
        % tensor de permeabilidade do elemento a esquerda
        
        Klef(1,1)=kmap(elem(lef,5),2);
        Klef(1,2)=kmap(elem(lef,5),3);
        Klef(2,1)=kmap(elem(lef,5),4);
        Klef(2,2)=kmap(elem(lef,5),5);
        
        % tensor de permeabilidade do elemento a direita
        
        Krel(1,1)=kmap(elem(rel,5),2);
        Krel(1,2)=kmap(elem(rel,5),3);
        Krel(2,1)=kmap(elem(rel,5),4);
        Krel(2,2)=kmap(elem(rel,5),5);
        %Determinação dos centróides dos elementos à direita e à esquerda.%
        C1 = centelem(inedge(iface,3),:); % baricentro do elemento a esquerda
        C2 = centelem(inedge(iface,4),:); % baricentro do elemento direito
        
        vd1 = coord(inedge(iface,2),:) - coord(inedge(iface,1),:);
        
        %Determinação das alturas dos centróides dos elementos à direita e à%
        %esquerda.                                                          %
        
        vd2 = C2 - coord(inedge(iface,1),:);     %Do início da aresta até o
        %centro da célula da direita.
        cd = cross(vd1,vd2);
        H2 = norm(cd)/norm(vd1); % altura a direita
        
        ve2 = C1 - coord(inedge(iface,1),:);
        
        ce = cross(vd1,ve2);
        H1 = norm(ce)/norm(vd1); % altura a esquerda
        
        Kn1 = (RotH(vd1)'*Klef*RotH(vd1))/norm(vd1)^2;
        Kn2 = (RotH(vd1)'*Krel*RotH(vd1))/norm(vd1)^2;
        
        % calculo das constantes nas faces internas
        Kde = ((Kn1*Kn2))/(Kn1*H2 + Kn2*H1);
        grav1=-(Klef*gravelem(lef,:)');
        gravface1= (H1/Kn1)*dot(((RotH(vd1)')'), grav1);
        grav2=-(Krel*gravelem(rel,:)');
        gravface2= (H2/Kn2)*dot(((-RotH(vd1)')'), grav2);
        %g(iface+size(bedge,1),1)=-Kde*(gravface1-gravface2);
        
         g(iface+size(bedge,1),1)=0.5*(dot(normals(iface+size(bedge,1),:)*Klef,gravelem(lef,:))+dot(normals(iface+size(bedge,1),:)*Krel,gravelem(rel,:)));
    end
        
    
    G(lef,1)=G(lef,1)-g(iface+size(bedge,1),1);
    G(rel,1)=G(rel,1)+g(iface+size(bedge,1),1);
    
end
end

function [RH]=RotH(vi)
%Função que retorna a rotação de um certo vetor em 90 graus,
%anti-horariamente, na primeira coluna da matriz R, e horariamente,
%na segunda. Restringe um vetor de 3 coordenadas a um de 2, considerando
%que a terceira coordenada é nula.
% vi2=zeros(2,1);
% if size(vi)~=[3 1]
%     vi=vi';
%     vi2(1)=vi(1);
%     vi2(2)=vi(2);
% end
RH=[0 1 0;-1 0 0;0 0 0]*vi';
end
