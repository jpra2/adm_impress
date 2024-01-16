%funçao que calcula os fluxos nas arestas internas
%equacoes 28 e 29 (heterogeneo) ou 15 e 16 (homogeneo)

function [flowrate, flowresult]=flowrateTPFA(p,Kde,Kn,Hesq,nflagface,mobility,gravresult,gravrate,pinterp)

global coord bedge inedge bcflag centelem


%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "bedgeamount"
bedgeamount = 1:bedgesize;

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
for ifacont=1:size(bedge,1);
    lef=bedge(ifacont,3);
    nor=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*nor);
    if bedge(ifacont,5)<200 % se os nós esteverem na fronteira de DIRICHLET
        c1=nflagface(ifacont,2);
        
        
        flowrate(ifacont)=mobility*A*(nor^2)*(c1-p(lef));
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        flowrate(ifacont)= nor*bcflag(r,2);
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3); %indice do elemento a direita da aresta i
    rel=inedge(iface,4); %indice do elemento a esquerda da aresta i
    
    %-------------------- calculo das vazões e velocidades ---------------%
    
    flowrate(iface+size(bedge,1))=mobility*Kde(iface)*(p(rel)-p(lef));
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
end

end