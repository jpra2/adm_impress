function [bc_matrix] = bcGLS()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global bedge coord bcflag
num_nodes = size(coord,1);
bc_matrix = sparse(num_nodes,8);
bnodes = unique(bedge(:,1:2));
for ii = bnodes'
    nodes_around = nsurnG(ii);
    start_node = nodes_around(1);
    end_node = nodes_around(end);
    edge1_nodes = [start_node ii];
    edge2_nodes = [ii end_node];
    tmp = ismember(bedge(:,1:2),edge1_nodes);
    edges1 = find(all(tmp,2));
    tmp = ismember(bedge(:,1:2),edge2_nodes);
    edges2 = find(all(tmp,2));
    flag1 = bedge(edges1,end);
    flag2 = bedge(edges2,end);    
    if flag1 > 200
        p1 = 2;        
    else
        p1 = 1;
    end
    if flag2 > 200
        p2 = 2;
    else
        p2 = 1;
    end        
    flagvalue1 = bcflag(bcflag(:,1) == flag1 ,2);
    flagvalue2 = bcflag(bcflag(:,1) == flag2 ,2);    
    bc_matrix(ii,1) = p1;
    bc_matrix(ii,2) = p2;
    bc_matrix(ii,3) = flagvalue1;
    bc_matrix(ii,4) = flagvalue2; 
    bc_matrix(ii,5) = flag1;
    bc_matrix(ii,6) = flag2; 
    bc_matrix(ii,7) = edges1;
    bc_matrix(ii,8) = edges2;     
   end
end

