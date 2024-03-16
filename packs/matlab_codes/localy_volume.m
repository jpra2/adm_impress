function local_volume = localy_volume(list,xt)

% the area of the element (triangle) is half the cross product of edge
% vectors
area=@(u,v)abs(u(1)*v(2)-u(2)*v(1))/2;

nlist = size(list,1);

local_volume = 0;

for i = 2:(nlist-1)
    
    v1 = xt(list(i),:)' - xt(list(1),:)';
    v2 = xt(list(i+1),:)' - xt(list(1),:)';
    
    a1 = area(v2,v1);
    
    local_volume = local_volume + a1;
    
end

end
    