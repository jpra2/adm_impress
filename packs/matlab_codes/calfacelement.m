function [auxface]=calfacelement(ielem,ifacont,aux11,aux12,bedge,inedge,auxface,posit1,posit2)

if aux11>size(bedge,1)
    auxlef=inedge(aux11-size(bedge,1),3);
    auxrel=inedge(aux11-size(bedge,1),4);
    auxvetor=[auxlef auxrel];
    a=auxvetor==ielem;
    r=find(a==0);
    
    auxface(ifacont,posit1)=auxvetor(r);
    
else
   auxlef=bedge(aux11,3);
   
   auxface(ifacont,posit1)=auxlef;
end

if aux12>size(bedge,1)
    auxlef=inedge(aux12-size(bedge,1),3);
    auxrel=inedge(aux12-size(bedge,1),4);
    auxvetor=[auxlef auxrel];
    a=auxvetor==ielem;
    r=find(a==0);
    
    auxface(ifacont,posit2)=auxvetor(r);
    
else
   auxlef=bedge(aux12,3);
   
   auxface(ifacont,posit2)=auxlef;
end

end