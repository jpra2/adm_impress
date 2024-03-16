function [ node ] = nsurnG(node )
global nsurn1 nsurn2 
point = nsurn2(node)+1;
length = nsurn2(node+1) - nsurn2(node);
node = nsurn1(point :(point+length-1));
end

