function [omega]=aux_omegacalculate(E,EE,netas,no)

global  esurn2 nsurn2

O=zeros(esurn2(no+1)-esurn2(no),3); % vetor de baricentro na vizinhança do nó "ni".
nec=size(O,1);
P=zeros(nsurn2(no+1)-nsurn2(no),3);


if size(P,1)==size(O,1) %Se for um nó interno.
    for k=1:nec,
        
        
        if (k==1)&&(size(P,1)==size(O,1))
            
            omega(k,1)=((netas(k,1)+netas(size(O,1),2))*(E(k,1))/(EE(k,1)+EE(size(O,1),2)))+ ((netas(k,2)+netas(k+1,1))*(E(k,2))/(EE(k,2)+EE(k+1,1)));
            
        elseif (k==nec)&&(size(P,1)==size(O,1))
            
            omega(k,1)=((netas(k,1)+netas(k-1,2))*(E(k,1))/(EE(k,1)+EE(k-1,2)))+ ((netas(k,2)+netas(1,1))*(E(k,2))/(EE(k,2)+EE(1,1)));
            
        else
            omega(k,1)=((netas(k,1)+netas(k-1,2))*(E(k,1))/(EE(k,1)+EE(k-1,2)))+ ((netas(k,2)+netas(k+1,1))*(E(k,2))/(EE(k,2)+EE(k+1,1)));
        end
        
        
        
    end
else %Se for um nó do contorno.
    for k=1:nec,

        if (k==1)&&(size(P,1)~=size(O,1))&& k~=nec
            
            omega(k,1)=((netas(k,1))*(E(k,1))/(EE(k,1)))+ ((netas(k,2)+netas(k+1,1))*(E(k,2))/(EE(k,2)+EE(k+1,1)));
                        
        elseif (k==1)&&(size(P,1)~=size(O,1))&& k==nec
            
            omega(k,1)=((netas(k,1))*(E(k,1))/(EE(k,1)))+ ((netas(k,2))*(E(k,2))/(EE(k,2)));
            
        elseif (k==nec)&&(size(P,1)~=size(O,1))
            
            omega(k,1)=((netas(k,1)+netas(k-1,2))*(E(k,1))/(EE(k,1)+EE(k-1,2)))+ (netas(k,2)*(E(k,2))/(EE(k,2)));
        else
            omega(k,1)=((netas(k,1)+netas(k-1,2))*(E(k,1))/(EE(k,1)+EE(k-1,2)))+ ((netas(k,2)+netas(k+1,1))*(E(k,2))/(EE(k,2)+EE(k+1,1)));
            
        end
        
    end
end


end
