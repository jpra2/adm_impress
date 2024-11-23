function [Sat_max_min]=Saturation_max_min(S_old)
global elem esurn1 esurn2
Sat_max_min=zeros(size(elem,1),2);
for ielem =1:size(elem,1)
    n1=elem(ielem,1);
    n2=elem(ielem,2);
    n3=elem(ielem,3);
    r=1;
    for j=[n1,n2,n3]
        list=esurn1(esurn2(j)+1:esurn2(j+1));
         
        for k=1:size(list,1)
            a(k,r)=S_old(list(k));
        end
        r=r+1;
    end
    Sat_max_min(ielem,1)=max(max([a]));
    Sat_max_min(ielem,2)=min(min([a]));
end