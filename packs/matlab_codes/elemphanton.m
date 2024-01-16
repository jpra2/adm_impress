function [r,z,face_bedge]=elemphanton(N)
global bedge esurn1 esurn2
for i=1:size(bedge,1)
    
    v1=bedge(i,1);
    v=esurn1(esurn2(v1)+1:esurn2(v1+1));
    if size(v)==1
        r(v1,1:4)=v;
    else
        k1=size(v,1);
        r(v1,1:k1)=v(1:k1,1)';
        vinvert=fliplr(r(v1,1:k1));
        r(v1,k1+1:2*k1)=vinvert(1, 1:k1);
    end
    z(i,1:4)=[bedge(i,1) bedge(i,2) bedge(i,3) bedge(i,3)];
    k=0;
    for ii=1:size(N,2)
       if N(v1,ii)~=0
        k=k+1;
       end
    end
    if k==2
       face_bedge(v1,:)= [N(v1,1:k) N(v1,1:k)];
    else
       face_bedge(v1,:)= [N(v1,1:k) N(v1,2:k-1)]; 
    end

end
end