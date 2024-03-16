function [Transmicont]=auxassemblematrixcontourLINEAR(auxiface,auxparameter,ielem,normcont,weightDMP,bedge)
m=size(bedge,1);

auxweight=0;

if auxiface >m % quando auxiface pertece a face 
    auxlef=weightDMP(auxiface-m,3);
    auxrel=weightDMP(auxiface-m,4);
    
    weight1=weightDMP(auxiface-m,1);
    weight2=weightDMP(auxiface-m,2);
    
    if auxlef==ielem
        auxweight= weight2;
    elseif auxrel==ielem
        auxweight=weight1;
    end
    
    Transmicont= auxweight*auxparameter*normcont;

    
end

end