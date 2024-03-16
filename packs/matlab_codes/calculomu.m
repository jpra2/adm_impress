function [mulef, murel,partialfluxlef2,partialfluxrel2]=calculomu(lef, rel,Fface,valuemin,p,pinterp,iface,gamma,parameter)
global coord inedge bedge

norma= norm(coord(inedge(iface,1),:)-coord(inedge(iface,2),:));

%% Elemento a esquerda

iatual=iface+size(bedge,1);
partialfluxlef2=gamma*norma*parameter(1,2,iatual)*(p(lef)-pinterp(iatual))+ norma*parameter(1,1,iatual)*(p(lef)-pinterp(parameter(1,3,iatual)));

%% Elemento a direita

partialfluxrel2=gamma*norma*parameter(2,2,iatual)*(p(rel)-pinterp(iatual))+ norma*parameter(2,1,iatual)*(p(rel)-pinterp(parameter(2,3,iatual)));

%% Calculo das constantes (eq. 13) do artigo Gao e Wu, 2013

 mulef=(abs(partialfluxrel2)+valuemin)/(abs(partialfluxlef2)+abs(partialfluxrel2)+2*valuemin);
 murel=(abs(partialfluxlef2)+valuemin)/(abs(partialfluxlef2)+abs(partialfluxrel2)+2*valuemin);
% mulef=(abs(partialfluxrel2))/(abs(partialfluxlef2)+abs(partialfluxrel2)+2*valuemin);
% murel=(abs(partialfluxlef2))/(abs(partialfluxlef2)+abs(partialfluxrel2)+2*valuemin);
end