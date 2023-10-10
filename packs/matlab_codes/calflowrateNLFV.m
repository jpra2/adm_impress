function [influx]=calflowrateNLFV(p, pinterp, parameter,nflag,kmap,Fface,valuemin,gamma,mobility)
global inedge coord bedge bcflag elem centelem
influx=zeros(size(inedge,1)+size(bedge,1),1);
for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    b= nflag(ifacont,1)>200;
    if b==1
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        influx(ifacont,1)= normcont*bcflag(r,2);
    else
        ifacelef=parameter(1,4,ifacont);
        if ifacelef >size(bedge,1)
            auxlef=inedge(ifacelef-size(bedge,1),3);
            auxrel=inedge(ifacelef-size(bedge,1),4);
            [weight1, weight2]=weightnlfv(kmap,ifacelef-size(bedge,1));
            
            if auxlef==lef
                elemavaliar=auxrel;
                auxweight= weight2;
            elseif auxrel==lef
                elemavaliar=auxlef;
                auxweight=weight1;
            end
            
            Transmicont= auxweight*parameter(1,5,ifacont)*normcont;
            
            influx(ifacont,1)=influx(ifacont,1) + mobility(ifacont)*Transmicont*(p(lef)-p(elemavaliar));
        else
            
            influx(ifacont,1)=influx(ifacont,1) + mobility(ifacont)*parameter(1,1,ifacont)*normcont*(p(lef)-pinterp(ifacelef));
        end
        
        
        influx(ifacont,1)=influx(ifacont,1)+ mobility(ifacont)*normcont*parameter(1,2,ifacont)*(p(lef)-pinterp(ifacont));
        
    end
    
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    %% calculo dos pesos
    
    [weightlef, weightrel]=weightnlfv(kmap,iface);
    %% Calculo dos fluxos parciais (eq. 11)
    ifactual=iface+size(bedge,1);
    a=parameter(1,1,ifactual);
    c=parameter(1,3,ifactual);
    
    if a==0 && parameter(1,4,ifactual)==0
        
        partialfluxlef=norma*(parameter(1,2,ifactual)*(p(lef)-pinterp(parameter(1,5,ifactual)))+ ...
            parameter(1,3,ifactual)*(p(lef)-pinterp(parameter(1,6,ifactual))));
        
    elseif c==0
        partialfluxlef=norma*(gamma*parameter(1,1,ifactual)*(p(lef)-pinterp(parameter(1,4,ifactual)))+ ...
            parameter(1,2,ifactual)*(p(lef)-pinterp(parameter(1,5,ifactual))));
        
    end
    % Fluxo parcial do elemento a direita
    a1=parameter(2,1,ifactual);
    c1=parameter(2,3,ifactual);
    
    if a1==0 && parameter(2,4,ifactual)==0
        
        partialfluxrel=norma*(parameter(2,2,ifactual)*(p(rel)-pinterp(parameter(2,5,ifactual)))+ ...
            parameter(2,3,ifactual)*(p(rel)-pinterp(parameter(2,6,ifactual))));
    elseif c1==0
        partialfluxrel=norma*(gamma*parameter(2,1,ifactual)*(p(rel)-pinterp(parameter(2,4,ifactual)))+ ...
            parameter(2,2,ifactual)*(p(rel)-pinterp(parameter(2,5,ifactual))));
    end
    
    %% Calculo das constantes da não linearidade (eq. 13) do artigo Gao e Wu, (2013)
    if valuemin==0 || partialfluxlef*partialfluxrel>0
        mulef=abs(partialfluxrel)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
        
        murel=abs(partialfluxlef)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
    else
        mulef=(abs(partialfluxrel)+valuemin)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
        
        murel=(abs(partialfluxlef)+valuemin)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
    end
    mulef1=(abs(partialfluxrel)+valuemin)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
    
    murel1=(abs(partialfluxlef)+valuemin)/(abs(partialfluxlef)+abs(partialfluxrel)+2*valuemin);
    %% implementação do beta
    betalef=mulef*(1-sign(partialfluxlef*partialfluxrel));
    
    %% calculo da transmissibilidade
    alfa=(1-gamma)*norma*(mulef1*parameter(1,1,iface+size(bedge,1))*weightrel + murel1*parameter(2,1,iface+size(bedge,1))*weightlef);
    Transmilef=alfa + weightrel*gamma*betalef*parameter(1,1,ifactual)*norma;
    
    influx(iface+ size(bedge,1),1)=influx(iface+ size(bedge,1),1)+ mobility(ifactual)*Transmilef*(p(lef)-p(rel));
    
    ifacelef=parameter(1,5,ifactual);
    if ifacelef >size(bedge,1)
        auxlef=inedge(ifacelef-size(bedge,1),3);
        auxrel=inedge(ifacelef-size(bedge,1),4);
        
        [weight1,weight2]=weightnlfv(kmap,ifacelef-size(bedge,1));
        
        c=find([auxlef auxrel]~=lef);
        if c==1
            elemavaliarlef= auxlef;
            auxweight= weight1;
        else
            elemavaliarlef= auxrel;
            auxweight= weight2;
        end
        
        auxTransmilef= auxweight*betalef*parameter(1,2,ifactual)*norma;
        influx(iface+ size(bedge,1),1)= influx(iface+ size(bedge,1),1)+ mobility(ifactual)*auxTransmilef*(p(lef)-p(elemavaliarlef));
        
    else
        
        auxTransmilef= betalef*parameter(1,2,ifactual)*norma;
        influx(iface+ size(bedge,1),1)= influx(iface+ size(bedge,1),1)+mobility(ifactual)*auxTransmilef*(p(lef)-pinterp(ifacelef));
        
    end
    ifacelef2=parameter(1,6,ifactual);
    if ifacelef2~=0
        if  ifacelef2 >size(bedge,1)
            auxlef=inedge(ifacelef2-size(bedge,1),3);
            auxrel=inedge(ifacelef2-size(bedge,1),4);
            
            [weight1,weight2]=weightnlfv(kmap,ifacelef2-size(bedge,1));
            
            t=find([auxlef auxrel]~=lef);
            if t==1
                elemavaliarlef= auxlef;
                auxweight= weight1;
            else
                elemavaliarlef= auxrel;
                auxweight= weight2;
            end
            
            auxTransmilef= auxweight*betalef*parameter(1,3,ifactual)*norma;
            influx(iface+ size(bedge,1),1)= influx(iface+ size(bedge,1),1)+ mobility(ifactual)*auxTransmilef*(p(lef)-p(elemavaliarlef));

        else
            
            auxTransmilef= betalef*parameter(1,3,ifactual)*norma;
            
            influx(iface+ size(bedge,1),1)= influx(iface+ size(bedge,1),1)+mobility(ifactual)*auxTransmilef*(p(lef)-pinterp(ifacelef));
        end
    end
    
end
end