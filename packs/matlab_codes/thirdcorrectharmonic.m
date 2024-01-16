function [yy,weightDMP,raioaux]=thirdcorrectharmonic(kmap)

global inedge bedge coord elem centelem 
yy=zeros(size(inedge,1)+size(bedge,1),3);
Klef=zeros(3,3);
Krel=zeros(3,3);
RR=[0 1 0; -1 0 0;0 0 0];

for ifacont=1:size(bedge,1)
    
    yy(ifacont,:)= 0.5*( coord(bedge(ifacont,1),:)+ coord(bedge(ifacont,2),:));    
end
i=1;
contador=1;
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    % comprimento da face ou area
    R=0.5*norm(vd1);
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
    
    % calculo do ponto de intersecao: XA e XB do artigo Zhang Kobaise
    
    xA=calculXaXb(coord(inedge(iface,1),:),coord(inedge(iface,2),:),centelem(lef,:),(Klef*RR*(vd1')/norm(vd1))');
    xB=calculXaXb(coord(inedge(iface,1),:),coord(inedge(iface,2),:),centelem(rel,:),(Krel*RR*(-vd1')/norm(vd1))');
    
    % calculo dos parametros em cada face interna
    w1=norm(Klef*RR*(vd1')/norm(vd1));
    w2=norm(Krel*RR*(-vd1')/norm(vd1));
    L1=norm(xA-centelem(lef,1:2));
    L2=norm(xB-centelem(rel,1:2)) ;
    a=w1/L1; b=w2/L2;
    
    % calculo das constantes normais em cada face interna
    Knlef=dot(RR*vd1',Klef*(RR*vd1')/norm(vd1)^2);
    
    Knrel=dot((RR*(-vd1')),Krel*(RR*(-vd1'))/norm(vd1)^2);
    % calculo dos pontos armonicos
    
    harmo=((hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)'+...
        hlef*hrel*(Klef'-Krel')*(RR*(vd1/norm(vd1))'))/(hrel*Knlef+hlef*Knrel))';
   
    % ponto medio da face
    media=0.5*(coord(inedge(iface,1),:)+coord(inedge(iface,2),:)); 
    %=====================================================================%
    % equação 37 Zhang e Kobaise
    raio=norm(harmo- media)/R;
    % R' is called as Raux
    Raux=3*R;
    if norm(harmo- media)>Raux % harmonic point lying outside, 
                
        raioaux(i,1)=raio;
        i=i+1;
        %% escolha do ponto harmonico mediante um processo de otimizacao     
        options = optimoptions(@fmincon,'Algorithm','sqp');
        % objetive function: equation 39
        fun =@(t)sqrt((a*(xA(1,1)-t(1))+b*(xB(1,1)-t(1)))^2+(a*(xA(1,2)-t(2))+b*(xB(1,2)-t(2)))^2); 
        % objetive function: equation A-4
        %fun=@(t)sqrt((harmo1(1,1)-t(1))^2+(harmo1(1,2)-t(1))^2);
        % initial value
        %t0 = media(1,1:2);
        t0=[0,0];
        % unitdisk: is restriction function.
        [yy(iface+size(bedge,1),1:2),fval] = fmincon(fun,t0,[],[],[],[],[],[],@(t) unitdisk(t,Raux,media),options);
        % negative value z, then satifies the restriction 40
        
        if abs(norm(yy(iface+size(bedge,1),1:2)- media(1,1:2))-Raux)>1e-5
          disp('...The restricction of the equation 40 was violate...!!!')
        end
        % calculo dos pesos para interpolar a variavel no ponto em questão
        weightlef_1= a/(a+b); weightrel_1= 1-weightlef_1; % Eq. (17)
        weightlef_2 = 0;weightrel_2 = 0;
        weightlef = weightlef_1+weightlef_2;
        weightrel = weightrel_1+weightrel_2;
        % serve para armazenar o numero de faces que precisa de corregir
        contador=contador+1;
    else
        raioaux=0;
        % escolha o original ponto harmonico
        yy(iface+size(bedge,1),:)=harmo;
        % calculo dos pesos para interpolar a variavel no ponto em questão
        weightlef_1= (hrel*Knlef)/(hrel*Knlef+ hlef*Knrel); weightrel_1= 1-weightlef_1; % Eq. (17)
        
        weightlef_2 = 0;weightrel_2 = 0;
        weightlef = weightlef_1+weightlef_2;
        weightrel = weightrel_1+weightrel_2;
    
        
    end
    %% Calculo dos pesos
    
    weightDMP(iface,1)=weightlef;
    weightDMP(iface,2)=weightrel;
    weightDMP(iface,3)=lef;
    weightDMP(iface,4)=rel;
    
   %======================================================================% 
end
name=(contador/size(inedge,1))*100;
sprintf('>> Percentage of faces corrected: %s ',num2str(name))
end

function [p]=calculXaXb(p1,p2,p3,p4)

A=[p1(1,2)-p2(1,2) p2(1,1)-p1(1,1);
    p3(1,2)-p4(1,2) p3(1,1)-p4(1,1)];

b=[p1(1,2)*p2(1,1)-p1(1,1)*p2(1,2);
   p3(1,2)*p4(1,1)-p3(1,1)*p4(1,2)];

p=(inv(A)*b)';
end

