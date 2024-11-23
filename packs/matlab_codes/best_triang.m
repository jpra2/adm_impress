function best_tri = best_triang(utiltri,xt,coord,node)

n = size(utiltri,1);
dist = zeros(1,n);
dist_min = zeros(n,1);
best_tri = zeros(1,3);

for i = 1:n
    
    % conectivity
    v1 = utiltri(i,1);
    v2 = utiltri(i,2);
    v3 = utiltri(i,3);
    
    % point coordinates
    coordv1 = xt(v1,:); %coordinates of point 1
    coordv2 = xt(v2,:); %coordenates of point 2
    coordv3 = xt(v3,:); %coordenates of point 3
    coordvn = coord(node,1:2); %coordinates of node
    
    % distances
    dist(1,i) = (norm(coordvn - coordv1) + norm(coordvn - coordv2) + ...
        norm(coordvn - coordv3));
end

dist_mini = sort(dist);
dist_min = dist_mini(1);

[a,b] = find(dist == dist_min);
best_tri(1,:) = utiltri(b,:);

end