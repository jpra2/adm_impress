function [M,I]=assemblematrixLPSPPS(p,pinterp,parameter,fonte)
global inedge coord bedge bcflag elem 
I=sparse(size(elem,1),1);
M=sparse(size(elem,1),size(elem,1));
valuemin=1e-16;
%% fonte
I=I+fonte;
%%
for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
  
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef)- normcont*bcflag(r,2);
    else  
        %% calculo da contribuição do contorno, veja Eq. 2.17 (resp. eq. 24) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)   
           alef=normcont*(parameter(1,1,ifacont)*pinterp(parameter(1,3,ifacont))+...
                parameter(1,2,ifacont)*pinterp(parameter(1,4,ifacont)));
            
            Alef=normcont*(parameter(1,1,ifacont)+parameter(1,2,ifacont));
            
            %% implementação da matriz global no contorno
            M(lef,lef)=M(lef,lef)+ Alef;
            I(lef,1)=I(lef,1)+alef;
    end
end
%% Montagem da matriz global

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    
    %% calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    % esquerda
    alef=parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
        parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual));
    % direita
    
    arel= parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
        parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual));
    %% calculo dos "mu", Eq. 2.8 (resp. eq. 18) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
     mulef=(abs(arel)+valuemin)/(abs(alef)+abs(arel)+2*valuemin);  
     murel=(abs(alef)+valuemin)/(abs(alef)+abs(arel)+2*valuemin);
    %% calculo da contribuição, Eq. 2.10 (resp. eq. 19) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    Bsigma=murel*arel-mulef*alef;
    
    Bmas=(abs(Bsigma)+Bsigma)/2;
    Bmenos=(abs(Bsigma)-Bsigma)/2;
    
    %% calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*(mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual))+ Bmas/(p(lef)+valuemin));
  %  ALR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual))+ Bmenoslef/(p(rel)+valuemin);
    
    ARR=norma*(murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual))+ Bmenos/(p(rel)+valuemin));
  %  ARL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual))+ Bmenosrel/(p(lef)+valuemin);
    %% implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ ALL;
    M(lef,rel)=M(lef,rel)- ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ ARR;
    M(rel,lef)=M(rel,lef)- ALL; 
    
    I(lef, 1)=I(lef,1)- (norma*Bmas*valuemin/(p(lef)+valuemin)-norma*Bmenos*valuemin/(p(rel)+valuemin));
    I(rel, 1)=I(rel,1)+ (norma*Bmas*valuemin/(p(lef)+valuemin)-norma*Bmenos*valuemin/(p(rel)+valuemin));
    
end
%% malha 23x23
% M(357,:)=0*M(357,:);
% M(357,357)=1;
% I(357)=1;
% M(173,:)=0*M(173,:);
% M(173,173)=1;
% I(173)=0;
%% malha 11x11
% M(83,:)=0*M(83,:);
% M(83,83)=1;
% I(83)=1;
% M(39,:)=0*M(39,:);
% M(39,39)=1;
% I(39)=0;
end