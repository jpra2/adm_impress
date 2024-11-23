function [ node ] = esurnG(node )
global esurn1 esurn2 
point = esurn2(node)+1;
length = esurn2(node+1) - esurn2(node);
node = esurn1(point :(point+length-1));



end