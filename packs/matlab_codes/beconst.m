function [BI,BJ] = beconst(bedge,coord,kmap,elem,xt,Pnode,vN,lambda_bedge)

area=@(u,v)abs(u(1)*v(2)-u(2)*v(1))/2;

Rot=-[cos(pi/2) -sin(pi/2); %matriz de rotaçao de pi/2 (horária)
      sin(pi/2)  cos(pi/2)]; 

BI = zeros(size(bedge,1),1);
BJ = BI;

for i=1:size(bedge,1)
    v1 = bedge(i,1);        %node 1
    v2 = bedge(i,2);        %node 2
    
    coordv1 = coord(v1,1:2)'; %coordinate of node 1
    coordv2 = coord(v2,1:2)'; %coordinate of node 2
    
    edge = coordv2 - coordv1; %edge e
    Ne = -Rot*edge;% normal to edge
    
    edgetype = bedge(i,5);
    
    Iele = bedge(i,3);%indice do elemento a direita da aresta i
     
    if (edgetype <= 200) % edge with Dirichlet condition
        
        KI = [kmap(elem(Iele,5),2) kmap(elem(Iele,5),3); ...
            kmap(elem(Iele,5),4) kmap(elem(Iele,5),5)];
        
        DI = lambda_bedge(i)*KI;
        
        %normals to sub-element
        r1 = xt(Iele,:)' - coordv1;
        Nr1 = Rot*r1;
        r2 = coordv2 - xt(Iele,:)';
        Nr2 = Rot*r2;

        %area of sub-element
        t = area(r1,edge);

        %coeficients
        BI(i) = -(((Nr1 + Nr2)'*(DI*Ne))/(2*t));
        BJ(i) = ((Pnode(v1)*Nr2 + Pnode(v2)*Nr1)'*(DI*Ne))/(2*t);
   else
        BI(i) = 0;
        BJ(i) = vN(i);
    end
    
end