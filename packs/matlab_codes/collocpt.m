function xt = collocpt(elem,coord,kmap)

% function that compute the collocation points inside the element
% This collocation points is not the barycentre, but it is set depending of
% the permeability tensor

nelem = size(elem,1);

Rot=-[cos(pi/2) -sin(pi/2) 0;         %clockwise rotation matrix of pi/2
    sin(pi/2)  cos(pi/2) 0; 0 0 0];

xt = zeros(nelem,3);

for i=1:nelem
    
    k = [kmap(elem(i,5),2) kmap(elem(i,5),3) 0; kmap(elem(i,5),4)...
        kmap(elem(i,5),5) 0; 0 0 0];
    
    % vertices of the triangle
    v1 = elem(i,1);
    v2 = elem(i,2);
    v3 = elem(i,3);

    % nodal coordinates
    coordv1 = coord(v1,:)';     %node 1
    coordv2 = coord(v2,:)';     %node 2
    coordv3 = coord(v3,:)';     %node 3

    % edge vectors
    edge1 = coordv3 - coordv2;           %face opposite to vertex 1
    edge2 = coordv1 - coordv3;           %face opposite to vertex 2
    edge3 = coordv2 - coordv1;           %face opposite to vertex 3
        
    % normal vectors pointing away from the element
    N1 = Rot*edge1;          %normal to edge1
    N2 = Rot*edge2;          %normal to edge2
    N3 = Rot*edge3;          %normal to edge3
           
    dN(1)=((k*N1)'*N1)^0.5; %lenght of vector n1 in metric induced by the
                            %permeability tensor
    dN(2)=((k*N2)'*N2)^0.5; %lenght of vector n2 in metric induced by the 
                            %permeability tensor
    dN(3)=((k*N3)'*N3)^0.5; %lenght of vector n3 in metric induced by the 
                            %permeability tensor
        
    s = sum(dN);            %denominator of weights
    
    lambda=(1/s)*dN;        %weights
    
    % collocation points
    xt(i,:) = lambda(1)*coordv1' + lambda(2)*coordv2' + lambda(3)*coordv3';
end

end
