function [w] = weights(A,esuel,esuel_coord)
global centelem

[m,n,p] = size(A);
% compute QR using Gram-Schmidt


for k=1:p
    for j = 1:n
        v = A(:,j,k);
        for i=1:j-1
            R(i,j,k) = Q(:,i,k)'*A(:,j,k);
            v = v - R(i,j,k)*Q(:,i,k);
        end
        R(j,j,k) = norm(v);
        Q(:,j,k) = v/R(j,j,k);
    end
  Q(:,3,:)=0;  
    beta= (R(1,2,k)*R(2,3,k) - R(1,3,k)*R(2,2,k))/(R(1,1,k)*R(2,2,k));
   % calculo dos pesos segum o livro Computational fluid Dynamic author Blazek
    for j=1:m
        alfa1=(esuel_coord(j,1,k)-centelem(k,1))/R(1,1,k)^2; % x
        alfa2=(esuel_coord(j,2,k)-centelem(k,2)-(R(1,2,k)/R(1,1,k))*(esuel_coord(j,1,k)-centelem(k,1)))/R(2,2,k)^2; % y
        %alfa3=(esuel_coord(j,3,k)-c(k,3)-(R(2,3,k)/R(2,2,k))*(esuel_coord(j,2,k)-c(k,2))+beta*(esuel_coord(j,1,k)-c(k,1)))/R(3,3,k)^2; % z
      
       % w(:,j,k)=[alfa1; alfa2 ; alfa3 ];
       x=alfa1 -(R(1,2,k)/R(1,1,k))*alfa2;% + beta*alfa3;
       y=alfa2 ;%-(R(2,3,k)/R(2,2,k))*alfa3;
       %z=alfa3;
       w(:,j,k)=[x; y; 0 ];
       
    end
        
    w(3,:,k)=0;    
end
Q(:,3,:)=0;

end
