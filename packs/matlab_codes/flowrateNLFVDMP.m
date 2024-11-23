function [flowrate,flowresult]=flowrateNLFVDMP(p, pinterp, parameter,nflag,kmap,gamma,weightDMP,mobility)
global inedge coord bedge bcflag centelem phasekey smethod

valuemin=1e-16;
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
for ifacont=1:size(bedge,1)
    
    %Define "mobonface" (for "bedge")
    %It is a One-phase flow. In this case, "mobility" is "1"
   
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    b= nflag(ifacont,1)>200;
    if b==1
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        flowrate(ifacont,1)= -normcont*bcflag(r,2);
    else
        ifacelef1= parameter(1,3,ifacont);
        ifacelef2= parameter(1,4,ifacont);
        flowrate(ifacont,1)= mobility*normcont*(parameter(1,1,ifacont)*...
            (p(lef)-pinterp(ifacelef1))+parameter(1,2,ifacont)*(p(lef)-pinterp(ifacelef2)));
        
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont,1);
end
x=1000;
y=1000;
for iface=1:size(inedge,1)
    
   
    
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    %% Calculo dos fluxos parciais
    ifactual=iface+size(bedge,1);
    % Fluxo parcial do elemento a esquerda
    
    F1=norma*(parameter(1,1,ifactual)*(p(lef)-pinterp(parameter(1,3,ifactual)))+ parameter(1,2,ifactual)*(p(lef)-pinterp(parameter(1,4,ifactual))));
    
    % Fluxo parcial do elemento a direita
    F2=norma*(parameter(2,1,ifactual)*(p(rel)-pinterp(parameter(2,3,ifactual)))+ parameter(2,2,ifactual)*(p(rel)-pinterp(parameter(2,4,ifactual))));
    if abs(F1)<1e-20
        F1=0;
    end
    if abs(F2)<1e-20
        F2=0;
    end
    %% quando F1*F2<=0
    if F2*F1<0 || F2*F1==0
        
        mu1=(abs(F2)+valuemin)/(abs(F1)+abs(F2)+2*valuemin);
        
        mu2=(abs(F1)+valuemin)/(abs(F1)+abs(F2)+2*valuemin);
        
        % Calculo das contribuições do elemento a esquerda
        
        flowrate(iface+ size(bedge,1),1)= (mu1*F1-mu2*F2);
        
        
    else %% F1*F2>0
        [F1b,]=calfluxopartial2(parameter(1,3,ifactual),parameter(1,4,ifactual),...
            parameter(1,1,ifactual), parameter(1,2,ifactual), gamma,lef,pinterp,ifactual,p,norma);
        % aveliando no elemento a direita
        
        [F2b,]=calfluxopartial2(parameter(2,3,ifactual),parameter(2,4,ifactual),...
            parameter(2,1,ifactual), parameter(2,2,ifactual), gamma,rel,pinterp,ifactual,p,norma);
        
        %% Calculo das constantes da não linearidade (eq. 13) do artigo Gao e Wu, (2013)
        % veja a REMARK 3.2 pag. 313 do mesmo artigo
        if F1b*F2b>0
            mulef=abs(F2b)/(abs(F1b)+abs(F2b)+2*valuemin);
            
            murel=abs(F1b)/(abs(F1b)+abs(F2b)+2*valuemin);
        else
            mulef=(abs(F2b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
            
            murel=(abs(F1b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        end
        
        
        mulef1=(abs(F2b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        
        murel1=(abs(F1b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        %
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
            x=0;
            
        end
        %------------------------------------------------------------%
        if ifacelef1==ifactual
            auxparameter1=parameter(1,1,ifactual);
            auxparameter11= parameter(1,2,ifactual);
            auxifacelef11=parameter(1,4,ifactual);
        elseif ifacelef2==ifactual
            
            auxparameter1=parameter(1,2,ifactual);
            auxparameter11= parameter(1,1,ifactual);
            auxifacelef11=parameter(1,3,ifactual);
        else
            auxparameter11= parameter(1,1,ifactual);
            auxifacelef11=parameter(1,3,ifactual);
            auxparameter22= parameter(1,2,ifactual);
            auxifacelef22=parameter(1,4,ifactual);
            auxparameter1=0;
            y=0;
            
        end
        %% implementação do beta
        betalef=mulef*(1-sign(F1b*F2b));
        %=================================================================%
        if F1b*F2b<0 || F1b*F2b==0
            %% Calculo das contribuições do elemento a esquerda          
            weightlef=weightDMP(ifactual-size(bedge,1),1);
            weightrel=weightDMP(ifactual-size(bedge,1),2);
            
            alfa=(1-gamma)*norma*(mulef1*auxparameter1*weightrel + murel1*auxparameter2*weightlef);
            
            flowrate(iface+ size(bedge,1),1)=alfa*(p(lef)-p(rel))+...
                betalef*norma*(gamma*auxparameter1*(p(lef)-pinterp(ifactual))+auxparameter11*(p(lef)-pinterp(auxifacelef11))) ;
            %flow(iface+ size(bedge,1),1)=mobility(ifactual)*alfa*(p(lef)-p(rel))+betalef*F1b ;
        elseif (F2b<0 && F1b<0 && x==0 && y==0) || (F2b>0 && F1b>0 && x==0 && y==0)
            
            flowrate(iface+ size(bedge,1),1)=norma*(auxparameter11*(p(lef)-pinterp(auxifacelef11))+auxparameter22*(p(lef)-pinterp(auxifacelef22))) ;
            
            %=================================================================%
            %% F1b*F2b>0
        else
                      
            weightlef=weightDMP(ifactual-size(bedge,1),1);
            weightrel=weightDMP(ifactual-size(bedge,1),2);
            
            % Calculo da transmisibilidade (eq. 18)
            alfa=(1-gamma)*norma*(mulef1*auxparameter1*weightrel + murel1*auxparameter2*weightlef);
            
            flowrate(iface+ size(bedge,1),1)=alfa*(p(lef)-p(rel));
            
        end
        
    end
    flowrate(iface+ size(bedge,1),1)=flowrate(iface+ size(bedge,1),1)*mobility;
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
end
%flowresult = flow_balance(flowrate);
end