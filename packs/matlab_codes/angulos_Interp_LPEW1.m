function [ fi2, fi1, theta2, theta1 ] = angulos_Interp_LPEW1( O, P, T, Qo, No)
%Retorna os ângulos fi e theta. Ex: fi1=[f1(cell1) fi1(cell2)...
%fi1(celln)], onde n=esurn2(ni+1)-esurn2(ni);
global esurn2
%Prealocação dos vetores.%

fi2=zeros(1,esurn2(No+1)-esurn2(No)); 
fi1=zeros(1,esurn2(No+1)-esurn2(No)); 
theta2=zeros(1,esurn2(No+1)-esurn2(No));
theta1=zeros(1,esurn2(No+1)-esurn2(No)); 



%Loop que percorre todos os elementos que concorrem no nó "ni".%

for k=1:size(fi2,2)
    
        
    %Determinação dos vetores necessários à obtenção dos cossenos:%
    
    % size(P,1)==size(O,1): quando o "ni" pertenece ao malha interna
    % k==size(fi2,2): pontero ao ultimo elemento
    v0=O(k,:)-Qo; 
    if (k==size(fi2,2))&&(size(P,1)==size(O,1)) 
        vfi2=T(1,:)-O(size(fi2,2),:);
        vth2=T(1,:)-Qo;
    else
       % quando "ni" esta na fronteira de Neumann ou Dirichlet
        vfi2=T(k+1,:)-O(k,:);
        vth2=T(k+1,:)-Qo;
    end
    vfi1=T(k,:)-O(k,:);
    vth1=T(k,:)-Qo;
    

    %Determinação dos ângulos:%
    fi2(k)=acos(dot(-v0,vfi2)/(norm(v0)*norm(vfi2)));
    fi1(k)=acos(dot(-v0,vfi1)/(norm(v0)*norm(vfi1)));
    theta2(k)=acos(dot(v0,vth2)/(norm(v0)*norm(vth2)));
    theta1(k)=acos(dot(v0,vth1)/(norm(v0)*norm(vth1)));
    

end

end

