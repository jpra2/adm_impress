function [ elem ] = reordnodelem( elem, centelem, coord )
%

for i=1:size(elem,1)
   v1=coord(elem(i,1),:)-centelem(i,:);
   if elem(i,4)==0
       n=3;
   else
       n=4;
   end
   for j=1:n
       v2=coord(elem(i,j),:)-centelem(i,:);
       bv=cross(v1,v2);
       b=bv(3);
       theta(j)=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
       if b<0
          theta(j)=(2*pi)-theta(j); 
       end
   end 
   ord=zeros(size(theta,2),1);
   theta1=theta;
   for k=1:size(theta,2)
       pos=find(theta1==min(theta1));
       e1(k)=elem(i,pos);
       theta1(pos)=2*pi;
   end
   elem(i,1:n)=e1;
   clear ord theta e1
end

end

