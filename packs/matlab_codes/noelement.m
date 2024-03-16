function [A]= noelement
global esurn1 esurn2 coord bedge
A=zeros(size(coord,1),4);

for j=1:size(coord,1)
    v=bedge(:,1)==j;
    a=max(v);
    if a==0
        vetor_elem=esurn1(esurn2(j)+1:esurn2(j+1));
        A(j,:)=vetor_elem';
    end
end

end