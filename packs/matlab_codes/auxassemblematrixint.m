function [flag1,flag2,flag3,auxweight,pressureface]=auxassemblematrixint(auxiface,ielem1,pinterp,ifactual,weightDMP,bedge)
nbedgetotal=size(bedge,1);
flag1=0;
flag2=0;
flag3=0;
pressureface=0;
auxweight=0;
if auxiface >nbedgetotal && auxiface==ifactual
    flag1=1;
    
elseif auxiface >nbedgetotal && auxiface~=ifactual
    
    auxlef=weightDMP(auxiface-nbedgetotal,3);
    auxrel=weightDMP(auxiface-nbedgetotal,4);
    
    weight1=weightDMP(auxiface-nbedgetotal,1);
    weight2=weightDMP(auxiface-nbedgetotal,2);
    
    if auxlef==ielem1
        auxweight= weight2;
    elseif auxrel==ielem1
        auxweight= weight1;
    end
    flag2=1;
else
    
    pressureface=pinterp(auxiface);
    flag3=1;
    
end

end