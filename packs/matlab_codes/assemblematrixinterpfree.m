function [M,I]=assemblematrixinterpfree(valuemin,p,pinterp,parameter,fonte)
global inedge coord bedge bcflag elem elemarea
I=zeros(size(elem,1),1);
M=zeros(size(elem,1),size(elem,1));
%% fonte
I=I+fonte.*elemarea;
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
        
        alef= normcont*(parameter(1,1,ifacont)*pinterp(parameter(1,3,ifacont))+...
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
    
    material1=elem(lef,5);
    material2=elem(rel,5);
    % material homogeneo
    if material1==material2
        %% calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        % esquerda
        alef=norma*(parameter(1,1,ifactual)*p(parameter(1,3,ifactual))+...
            parameter(1,2,ifactual)*p(parameter(1,4,ifactual)));
        % direita
        
        arel= norma*(parameter(2,1,ifactual)*p(parameter(2,3,ifactual))+...
            parameter(2,2,ifactual)*p(parameter(2,4,ifactual)));
        %% calculo dos "mu", Eq. 2.8 (resp. eq. 18) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        if alef==0 && arel==0
            mulef= 0.5;
            murel=1-mulef;
        else
            mulef=abs(arel)/(abs(alef)+abs(arel));
            %mulef=arel/(alef+arel);
            murel=1-mulef;
        end
        %% calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        ALL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
        ALR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
        
        ARR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
        ARL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
        %% implementação da matriz global
        % contribuição da transmisibilidade no elemento esquerda
        M(lef,lef)=M(lef,lef)+ ALL;
        M(lef,rel)=M(lef,rel)- ALR;
        % contribuição da transmisibilidade no elemento direita
        M(rel,rel)=M(rel,rel)+ ARR;
        M(rel,lef)=M(rel,lef)- ARL;
    else
        %% Calculo dos fluxos parciais
    ifactual=iface+size(bedge,1);
    % aveliando no elemento a esquerda
    
    [partialfluxlef,parameterlef]=calfluxopartial2(parameter(1,3,ifactual),parameter(1,4,ifactual),...
        parameter(1,1,ifactual), parameter(1,2,ifactual), gamma,lef,pinterp,ifactual,p,norma);
    % aveliando no elemento a direita
    
    [partialfluxrel,parameterel]=calfluxopartial2(parameter(2,3,ifactual),parameter(2,4,ifactual),...
        parameter(2,1,ifactual), parameter(2,2,ifactual), gamma,rel,pinterp,ifactual,p,norma);
    
    %% Calculo das constantes da não linearidade (eq. 13) do artigo Gao e Wu, (2013)
    % veja a REMARK 3.2 pag. 313 do mesmo artigo
    if partialfluxlef*partialfluxrel>0
        mulef=abs(partialfluxrel)/(abs(partialfluxlef)+abs(partialfluxrel));
        
        murel=abs(partialfluxlef)/(abs(partialfluxlef)+abs(partialfluxrel));
    else
        mulef=(abs(partialfluxrel)+valuemin)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
        
        murel=(abs(partialfluxlef)+valuemin)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
    end
    %% implementação do beta
    betalef=mulef*(1-sign(partialfluxlef*partialfluxrel));
    betarel=murel*(1-sign(partialfluxlef*partialfluxrel));
    
    mulef1=(abs(partialfluxrel)+valuemin)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
    
    murel1=(abs(partialfluxlef)+valuemin)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
    
    if partialfluxrel*partialfluxlef<0 || partialfluxrel*partialfluxlef==0 || (abs(partialfluxrel)<1e-10 && abs(partialfluxlef)<1e-10)
        %% Calculo das contribuições do elemento a esquerda
        [weightlef,weightrel]=weightnlfvDMP(kmap,ifactual-size(bedge,1));
        ifacelef1=parameter(1,3,ifactual);
        ifacelef2=parameter(1,4,ifactual);
        ifacerel1= parameter(2,3,ifactual);
        ifacerel2= parameter(2,4,ifactual);
        if ifacerel1==ifactual
            auxparameter2=parameter(2,1,ifactual);
        elseif ifacerel2==ifactual
            
            auxparameter2=parameter(2,2,ifactual);
        else
            auxparameter2=0;
        end
        if ifacelef1==ifactual
            auxparameter1=parameter(1,1,ifactual);
        elseif ifacelef2==ifactual
            
            auxparameter1=parameter(1,2,ifactual);
        else
            auxparameter1=0;
        end
        
        alfa=(1-gamma)*norma*(mulef*auxparameter1*weightrel + murel*auxparameter2*weightlef);
        
        %% contribuição da transmisibilidade no elemento esquerda
        M(lef,lef)=M(lef,lef)+ alfa;
        M(lef,rel)=M(lef,rel)- alfa;
        
        M(rel,rel)=M(rel,rel)+ alfa;
        M(rel,lef)=M(rel,lef)- alfa;
        
        
        [M,I]=auxassemblematrixint(ifacelef1,M,I,parameter(1,1,ifactual),auxparameter2,...
            lef,rel,mulef1,murel1,weightlef,weightrel,pinterp,kmap,norma,betalef,ifactual,gamma);
        
        
        [M,I]=auxassemblematrixint(ifacelef2,M,I,parameter(1,2,ifactual),auxparameter2,...
            lef,rel,mulef1,murel1,weightlef,weightrel,pinterp,kmap,norma,betalef,ifactual,gamma);
        
        %% contribuição da transmisibilidade no elemento direita
        
        [M,I]=auxassemblematrixint(ifacerel1,M,I,parameter(2,1,ifactual),auxparameter1,...
            rel,lef,murel1,mulef1,weightrel,weightlef,pinterp,kmap,norma,betarel,ifactual,gamma);
        
        
        [M,I]=auxassemblematrixint(ifacerel2,M,I,parameter(2,2,ifactual),auxparameter1,...
            rel,lef,murel1,mulef1,weightrel,weightlef,pinterp,kmap,norma,betarel,ifactual,gamma);
    else
        ifacelef1=parameter(1,3,ifactual);
        ifacelef2=parameter(1,4,ifactual);
        ifacerel1= parameter(2,3,ifactual);
        ifacerel2= parameter(2,4,ifactual);
        if ifacerel1==ifactual
            auxparameter2=parameter(2,1,ifactual);
        elseif ifacerel2==ifactual
            
            auxparameter2=parameter(2,2,ifactual);
        else
            auxparameter2=0;
            y=0;
        end
        if ifacelef1==ifactual
            auxparameter1=parameter(1,1,ifactual);
        elseif ifacelef2==ifactual
            
            auxparameter1=parameter(1,2,ifactual);
        else
            auxparameter1=0;
            x=0;
        end
        
        %if x==0 && y==0
            
        %    disp('resolver aqui')
        %else
            [weightlef,weightrel]=weightnlfvDMP(kmap,ifactual-size(bedge,1));
            % Calculo da transmisibilidade (eq. 18)
            alfa=(1-gamma)*norma*(mulef1*auxparameter1*weightrel + murel1*auxparameter2*weightlef);
            if abs(alfa)<1e-10
                disp('matriz singular')
                iface
                alfaiface=alfa
                fluxpartialef=partialfluxlef
                fluxpartialrel=partialfluxrel
            end
            % contribuição da transmisibilidade no elemento esquerda
            M(lef,lef)=M(lef,lef)+ alfa;
            M(lef,rel)=M(lef,rel)- alfa;
            % contribuição da transmisibilidade no elemento direita
            M(rel,rel)=M(rel,rel)+ alfa;
            M(rel,lef)=M(rel,lef)- alfa;
        %end
    end
        
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
end