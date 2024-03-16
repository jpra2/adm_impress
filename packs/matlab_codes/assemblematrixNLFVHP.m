function [M,I]=assemblematrixNLFVHP(pinterp,parameter,fonte,wells,...
    Hesq,Kn,Kt,nflag,gravrate,gravresult,nflagface)
global inedge coord bedge bcflag elem elemarea centelem gravitational strategy
I=sparse(size(elem,1),1);
M=sparse(size(elem,1),size(elem,1));

%% fonte
I=I+fonte;
%%
m=0;
for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef)- normcont*bcflag(r,2);
    else
        
        %do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        if strcmp(gravitational,'yes')
            if strcmp(strategy,'starnoni')
                m=gravrate(ifacont);
            elseif strcmp(strategy,'inhouse')
                m=0;
            end
        end
        % calculo da contribuição do contorno, veja Eq. 2.17 (resp. eq. 24)
        alef= normcont*(parameter(1,1,ifacont)*pinterp(parameter(1,3,ifacont))+...
            parameter(1,2,ifacont)*pinterp(parameter(1,4,ifacont)));
        
        Alef=normcont*(parameter(1,1,ifacont)+parameter(1,2,ifacont));
        % implementação da matriz global no contorno
        M(lef,lef)=M(lef,lef)+ Alef;
        I(lef)=I(lef)+alef+m;
    end
    
end
%% Montagem da matriz global
%coef=max(calnormface)^2;
coef=1e-16;
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma= sqrt(vd1(1,1)^2+vd1(1,2)^2);
    ifactual=iface+size(bedge,1);
    
    % esquerda
    
    alef=parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
        parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual));
    
    % direita
    
    arel= parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
        parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual));
    
    mulef=(abs(arel)+coef)/(abs(alef)+abs(arel)+2*coef);
    murel=(abs(alef)+coef)/(abs(alef)+abs(arel)+2*coef);
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and
    %Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    
    ARR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ ALL;
    M(lef,rel)=M(lef,rel)- ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ ARR;
    M(rel,lef)=M(rel,lef)- ALL;
    
    if strcmp(gravitational,'yes')
        if strcmp(strategy,'starnoni')
            m=gravrate(size(bedge,1)+iface,1);
        elseif strcmp(strategy,'inhouse')
            m=0;
        end
        I(lef)=I(lef)+m;
        I(rel)=I(rel)-m;
    end
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
% for ielem=1:size(elem,1)
%     if elem(ielem,5)==4
%         M(ielem,:)=0*M(ielem,:);
%         M(ielem,ielem)=1;
%         I(ielem)=300;
%     elseif elem(ielem,5)==5
%          M(ielem,:)=0*M(ielem,:);
%         M(ielem,ielem)=1;
%         I(ielem)=300;
%     elseif elem(ielem,5)==6
%          M(ielem,:)=0*M(ielem,:);
%         M(ielem,ielem)=1;
%         I(ielem)=300;
%     end
% end
% adequação da matriz nos poços produtores
% if max(max(wells))~=0
%     for iw = 1:size(wells,1)
%         if wells(iw,2)==2 %produtor
%             M(wells(iw,1),:)=0*M(wells(iw,1),:);
%             M(wells(iw,1),wells(iw,1))=1;
%             I(wells(iw,1))=0;
%         end
%     end
% end
end