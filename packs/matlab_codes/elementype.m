function [k]=elementype(ielem)
global elem

k=0;

for i=1:size(elem,2)-1
    
    if elem(ielem,i)~=0
        k=k+1;
    end
end
end