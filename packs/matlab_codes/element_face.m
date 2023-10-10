function [F]=element_face
global bedge inedge elem coord
for ii=1:size(elem,1)
    i=1;
    n=elem(ii,4);
    if n~=0
        list=elem(ii,1:4);
    else
        list=elem(ii,1:3);
    end
    for jj=1:length(list)
        if jj<length(list)
            ibedge=find(((bedge(:,1)==list(jj+1) & bedge(:,2)==list(jj))|(bedge(:,1)==list(jj) & bedge(:,2)==list(jj+1))));
            
            if ibedge~=0
                F(ii,i)=ibedge; % flag da face 
                i=i+1;
            else
                iedge=find((inedge(:,1)==list(jj+1) & inedge(:,2)==list(jj))|(inedge(:,1)==list(jj) & inedge(:,2)==list(jj+1)));
                F(ii,i)=iedge+size(bedge,1);
                i=i+1;
            end
        else
            ibedge=find(((bedge(:,1)==list(jj) & bedge(:,2)==list(1))|(bedge(:,1)==list(1) & bedge(:,2)==list(jj))));
            
            if ibedge~=0
                F(ii,i)=ibedge;
                i=i+1;
            else
                iedge=find((inedge(:,1)==list(jj) & inedge(:,2)==list(1))|(inedge(:,1)==list(1) & inedge(:,2)==list(jj)));
                F(ii,i)=iedge+size(bedge,1);
                

                i=i+1;
            end
        end
    end
    
    
end
end