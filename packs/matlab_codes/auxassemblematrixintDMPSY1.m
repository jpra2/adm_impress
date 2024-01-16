function [flag1,flag2,auxweight,pressureface]=auxassemblematrixintDMPSY1(auxiface,ielem,pinterp,weightDMP,bedge)

nbedgetotal=size(bedge,1);
flag1=0;
flag2=0;
pressureface=0;
if auxiface >nbedgetotal
    % calcula os elementos que compartilham a face interior "auxiface"
    
    auxlef=weightDMP(auxiface-nbedgetotal,3);
    auxrel=weightDMP(auxiface-nbedgetotal,4);
    % calcula os pesos correpondente a os elementos "auxlef" e "auxrel"
    weight1=weightDMP(auxiface-nbedgetotal,1);
    weight2=weightDMP(auxiface-nbedgetotal,2);
    
    if auxlef==ielem
        auxweight= weight2;
    elseif auxrel==ielem
        auxweight= weight1;
    end
    flag1=1;
else
    auxweight=0;
    flag2=1;
    pressureface=pinterp(auxiface);
    
end

