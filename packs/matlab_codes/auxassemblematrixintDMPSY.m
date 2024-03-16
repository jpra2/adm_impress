function [M,I]=auxassemblematrixintDMPSY(auxiface,M,I,parameter1,parameter2,ielem1,ielem2,mu1,mu2,weight1,weight2,...
                                       pinterp,kmap,norma,beta,ifactual,gamma)
global bedge inedge 
if auxiface >size(bedge,1) && auxiface==ifactual
       
    Transmi1=weight2*gamma*beta*parameter1*norma;
    %% contribuição da transmisibilidade no elemento esquerda
    M(ielem1,ielem1)=M(ielem1,ielem1)+ Transmi1;
    M(ielem1,ielem2)=M(ielem1,ielem2)- Transmi1;
    
elseif auxiface >size(bedge,1) && auxiface~=ifactual
    
    auxlef=inedge(auxiface-size(bedge,1),3);
    auxrel=inedge(auxiface-size(bedge,1),4);
    
    [weight1,weight2]=weightnlfvDMP(kmap,auxiface-size(bedge,1));
    
    if auxlef==ielem1
        elemavaliarlef= auxrel;
        auxweight= weight2;
    elseif auxrel==ielem1
        elemavaliarlef= auxlef;
        auxweight= weight1;
    end
    
    auxTransmilef= auxweight*beta*parameter1*norma;
    
    M(ielem1,ielem1)=M(ielem1,ielem1)+ auxTransmilef;
    
    M(ielem1,elemavaliarlef)=  M(ielem1,elemavaliarlef)- auxTransmilef;  
    
else
    
    auxTransmilef= beta*parameter1*norma;
    
    M(ielem1,ielem1)=M(ielem1,ielem1)+ auxTransmilef;
    
    I(ielem1)=I(ielem1)+ auxTransmilef*pinterp(auxiface);
end

