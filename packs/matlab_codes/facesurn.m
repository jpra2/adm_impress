function [N]=facesurn
global nsurn1 nsurn2 inedge coord bedge
for no=1:size(coord,1)
    m=1;
    vetor1=nsurn1(nsurn2(no)+1:nsurn2(no+1));
    V=vetor1';
    for j= [V]
        ibedge=find(((bedge(:,1)==no & bedge(:,2)==j)|(bedge(:,1)==j & bedge(:,2)==no)));
        if ibedge~=0
            
            N(no,m)=ibedge+size(inedge,1);
            m=m+1;
        else
            
            iedge=find(((inedge(:,1)==j & inedge(:,2)==no)|(inedge(:,1)==no & inedge(:,2)==j)));
            
            N(no,m)=iedge;
            m=m+1;
        end
    end
    
    
end