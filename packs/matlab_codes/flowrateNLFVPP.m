function [flowrate,flowresult,coercividade]=flowrateNLFVPP(p, pinterp, parameter,mobility)
global inedge coord bedge bcflag centelem phasekey smethod
valuemin=1e-16;
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "bedgeamount"
bedgeamount = 1:bedgesize;

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);

denomi=0;
numera=0;
for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    %======================================================================
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    C1=centelem(lef,:); % baricentro do elemento a esquerda
    vd1=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    %Determinação das alturas dos centróides dos elementos à direita e à%
    %esquerda.                                                          %
    
    ve2=C1-coord(bedge(ifacont,1),:);
    
    ce=cross(vd1,ve2);
    H1=norm(ce)/norm(vd1); % altura a esquerda
    
    %======================================================================
    
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        flowrate(ifacont,1)= normcont*bcflag(r,2);
    else
              
        
       flowrate(ifacont,1)=normcont*(parameter(1,1,ifacont)*(p(lef)-pinterp(parameter(1,3,ifacont)))+...
            parameter(1,2,ifacont)*(p(lef)-pinterp(parameter(1,4,ifacont))));
        % teste de coercividade
        
        
        denomi=denomi+(normcont/H1)*(p(lef)-pinterp(parameter(1,3,ifacont)))^2;
        
        numera=numera+flowrate(ifacont,1)*(p(lef)-pinterp(parameter(1,3,ifacont)));
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont,1);
    
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    %======================================================================
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    C1=centelem(inedge(iface,3),:); % baricentro do elemento a esquerda
    C2=centelem(inedge(iface,4),:); % baricentro do elemento direito
    %Determinação das alturas dos centróides dos elementos à direita e à%
    %esquerda.                                                          %
    
    vd2=C2-coord(inedge(iface,1),:);     %Do início da aresta até o
    %centro da célula da direita.
    cd=cross(vd1,vd2);
    H2=norm(cd)/norm(vd1); % altura a direita
    
    ve2=C1-coord(inedge(iface,1),:);
    
    ce=cross(vd1,ve2);
    H1=norm(ce)/norm(vd1); % altura a esquerda
    
    %======================================================================
    %% calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    % esquerda
    alef=(parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
        parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual)));
    % direita
    
    arel=(parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
        parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual)));
     %% calculo dos "mu", Eq. 2.8 (resp. eq. 18) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    
    mulef=(abs(arel)+valuemin)/(abs(alef)+abs(arel)+2*valuemin);
    murel=(abs(alef)+valuemin)/(abs(alef)+abs(arel)+2*valuemin);
    %% calculo da contribuição, Eq. 2.10 (resp. eq. 19) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    Bsigma=murel*arel-mulef*alef;
    
    Bmas=(abs(Bsigma)+Bsigma)/2;
    Bmenos=(abs(Bsigma)-Bsigma)/2;

    %% calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*(mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual))+ Bmas/(p(lef)+valuemin));
    ARR=norma*(murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual))+ Bmenos/(p(rel)+valuemin));
    
    flowrate(iface+size(bedge,1),1)=ALL*p(lef)-ARR*p(rel)+(Bmas*valuemin/(p(lef)+valuemin))-(Bmenos*valuemin/(p(rel)+valuemin));
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(iface+size(bedge,1),1);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(iface+size(bedge,1),1);
    
    % teste de coercividade
    
    denomi=denomi+(norma/(H1+H2))*(p(lef)-p(rel))^2;
    
    numera=numera+flowrate(iface+size(bedge,1),1)*(p(lef)-p(rel));
    
end
coercividade=numera/denomi;
end