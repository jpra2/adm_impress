function [g]=gravnode(N,kmap,No,gravelem)

global inedge bedge normals elem coord

faces=N(No,:);
R=[0 1 ; -1 0 ];
Klef=zeros(2,2);
Krel=zeros(2,2);
G=0;
j=1;
for i=faces
    
    
    if i~=0
        
        if i<=size(bedge,1)
            T=(coord(bedge(i,1),1:2)+coord(bedge(i,2),1:2))*0.5; 
            Normal=R*(coord(No,1:2)-T)';
            ielem=bedge(i,3);
            Klef(1,1)=kmap(elem(ielem,5),2);
            Klef(1,2)=kmap(elem(ielem,5),3);
            Klef(2,1)=kmap(elem(ielem,5),4);
            Klef(2,2)=kmap(elem(ielem,5),5);
            
            g(j)=dot(Normal'*Klef,(gravelem(ielem,1:2)));
            % G=G-g;
        else
            no1= inedge(i-size(bedge,1),1);
            no2= inedge(i-size(bedge,1),2);
            T=(coord(no1,1:2)+coord(no2,1:2))*0.5; 
            Normal=R*(coord(No,1:2)-T)';
            ielem1=inedge(i-size(bedge,1),3);
            ielem2=inedge(i-size(bedge,1),4);
            Klef(1,1)=kmap(elem(ielem1,5),2);
            Klef(1,2)=kmap(elem(ielem1,5),3);
            Klef(2,1)=kmap(elem(ielem1,5),4);
            Klef(2,2)=kmap(elem(ielem1,5),5); 
            
            Krel(1,1)=kmap(elem(ielem2,5),2);
            Krel(1,2)=kmap(elem(ielem2,5),3);
            Krel(2,1)=kmap(elem(ielem2,5),4);
            Krel(2,2)=kmap(elem(ielem2,5),5);
            
          g(j)=0.5*(dot(Normal'*Klef,gravelem(ielem1,1:2))+dot(Normal'*Krel,gravelem(ielem2,1:2)));  
         % G=G-g;
        end
        j=j+1;
    end
end