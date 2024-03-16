function [weight]=weightlfv(parameter)
global inedge bedge 

R=[0 1 0; -1 0 0;0 0 0];
for iface=1:size(inedge,1);
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    ifactual=iface+size(bedge,1);
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=parameter(1,1,ifactual)+parameter(1,2,ifactual);
    ARR=parameter(2,1,ifactual)+parameter(2,2,ifactual);

%% Calculo dos pesos
weight(iface,1)=ALL/(ALL+ARR);
weight(iface,2)=ARR/(ALL+ARR);
weight(iface,3)=lef;
weight(iface,4)=rel;
end
end