function [T]=faces_no(bedge, inedge,No,m,n)

T=zeros(2,1);
     i=1;
     for jj=[m,n]
      
            ibedge=find(((bedge(:,1)==No & bedge(:,2)==jj)|(bedge(:,1)==jj & bedge(:,2)==No)));
            
            %b=logical(((bedge(:,1)==No & bedge(:,2)==jj)|(bedge(:,1)==jj & bedge(:,2)==No)));
            %key=1:length(b);
            %ibedge=key(b==1);
            if ibedge~=0
                T(i,1)=ibedge + size(inedge,1);
                i=i+1;
            else
                iedge=find((inedge(:,1)==No & inedge(:,2)==jj)|(inedge(:,1)==jj & inedge(:,2)==No));
                %b=logical((inedge(:,1)==No & inedge(:,2)==jj)|(inedge(:,1)==jj & inedge(:,2)==No));
                %key=1:length(b);
                %iedge=key(b==1);
                T(i,1)=iedge;
                i=i+1;
            end
     end
end