function [flowrate,flowresult]=flowratelfvLPEW(parameter,weightDMP,mobility,pinterp,p)
global inedge coord bedge bcflag centelem
% incialização das matrizes

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
    
for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        flowrate(ifacont,1)=normcont*bcflag(r,2);
    else
        nolef1=parameter(1,3,ifacont);
        nolef2=parameter(1,4,ifacont);
        
        flowrate(ifacont,1)= mobility*normcont*((parameter(1,1,ifacont)+parameter(1,2,ifacont))*p(lef)-...
            parameter(1,1,ifacont)*pinterp(nolef1)-parameter(1,2,ifacont)*pinterp(nolef2));
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
    
end
% Montagem da matriz global

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1); % F
    murel=weightDMP(ifactual-size(bedge,1),2); % B
    %[mulef, murel]=weightnlfvDMP(kmap,iface);
    % os nós que conforman os pontos de interpolação no elemento a esquerda
    auxnolef1=parameter(1,3,ifactual);
    auxnolef2=parameter(1,4,ifactual);
    % os nós que conforman os pontos de interpolação no elemento a direita
    auxnorel1=parameter(2,3,ifactual);
    auxnorel2=parameter(2,4,ifactual);
    % calculo dos fluxo parcial a esquerda
    fluxesq=norma*((parameter(1,1,ifactual)+parameter(1,2,ifactual))*p(lef)-...
                   parameter(1,1,ifactual)*pinterp(auxnolef1)-parameter(1,2,ifactual)*pinterp(auxnolef2));
    % calculo dos fluxo parcial a direita        
    fluxdireit=norma*((parameter(2,1,ifactual)+parameter(2,2,ifactual))*p(rel)-...
                    parameter(2,1,ifactual)*pinterp(auxnorel1)-parameter(2,2,ifactual)*pinterp(auxnorel2));
    % calculo do fluxo unico na face
    flowrate(iface+size(bedge,1),1)=mobility*(murel*fluxesq-mulef*fluxdireit);
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);  
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
    
end

end