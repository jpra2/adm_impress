function [elemento]=searchelement(a,b,c,d)

global elem
for ielem=1:size(elem,1)
    
    
    j=1;
    p=[0,0,0,0];
%     for t=elem(ielem,1:4)
%     
%         if t==a || t==b ||t==c || t==d
%             
%             x=t;
%             
%         else
%             x=0;
%         end
%         if x==0
%             break
%             
%         else
%             p(j)=x;
%             j=j+1;
%             elemento=ielem;
%         end
%                
%     end
%     
    
     for t=elem(ielem,1:3)
    
        if t==a || t==b ||t==c
            
            x=t;
            
        else
            x=0;
        end
        if x==0
            break
            
        else
            p(j)=x;
            j=j+1;
            elemento=ielem;
        end
               
    end
    
    if max(p(1,1))~=0 & max(p(1,2))~=0 & max(p(1,3))~=0 & max(p(1,4))~=0
        break
    end
     
end

end