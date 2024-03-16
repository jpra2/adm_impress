function wells = mod_wells(center,tetha)
j=2;
wells(1,1)=1301;
 wells(1,2)=1;
% wells(2,1)=729;
% wells(2,2)=2;
% wells(3,1)= 1286;
% wells(3,2)=2;
for i=1:size(center,1)
    
    a=[-0.3*sin(tetha),-0.3*cos(tetha), 0]; % poço produtor 2
    
    b= [0.3*sin(tetha),-0.3*cos(tetha) , 0]; %poço produtor 1
    
    if norm(b(1,:)-center(i,:))<0.01 ||norm(b(1,:)-center(i,:))==0  
        wells(j,1)=i;
        wells(j,2)=2;
        j=j+1;
    elseif norm(a(1,:)-center(i,:))<0.01 || norm(a(1,:)-center(i,:))==0
        wells(j,1)=i;
        wells(j,2)=2;
        j=j+1;
    end
    
end

end