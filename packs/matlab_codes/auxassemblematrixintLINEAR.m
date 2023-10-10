function [auxweight,auxelemval]=auxassemblematrixintLINEAR(auxiface,ielem,weightDMP)

global bedge
nbedgetotal=size(bedge,1);
    % calcula os elementos que compartilham a face interior "auxiface"
    
    auxlef=weightDMP(auxiface-nbedgetotal,3);
    auxrel=weightDMP(auxiface-nbedgetotal,4);
    % calcula os pesos correpondente a os elementos "auxlef" e "auxrel"
    weight1=weightDMP(auxiface-nbedgetotal,1);
    weight2=weightDMP(auxiface-nbedgetotal,2);
    
    if auxlef==ielem
        auxweight= weight2;
        auxelemval=auxrel;
    elseif auxrel==ielem
        auxweight= weight1;
        auxelemval=auxlef;
    end
end