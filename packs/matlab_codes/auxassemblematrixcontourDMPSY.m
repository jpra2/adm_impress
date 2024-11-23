function [Transmicont,elemavaliar,pressureface,contorno]=...
    auxassemblematrixcontourDMPSY(auxiface,auxparameter,ielem,pinterp,normcont,weightDMP,bedge)
m=size(bedge,1);

if auxiface >m % quando auxiface pertece a face 
    auxlef=weightDMP(auxiface-m,3);
    auxrel=weightDMP(auxiface-m,4);
    
    weight1=weightDMP(auxiface-m,1);
    weight2=weightDMP(auxiface-m,2);
    
    if auxlef==ielem
        elemavaliar=auxrel;
        auxweight= weight2;
    elseif auxrel==ielem
        elemavaliar=auxlef;
        auxweight=weight1;
    end
    
    Transmicont= auxweight*auxparameter*normcont;
    
    contorno=0;
    pressureface=0;
    
else
    
    Transmicont= auxparameter*normcont;
    
    contorno=1;
    
    pressureface=pinterp(auxiface);
    
    elemavaliar=ielem;

end

end