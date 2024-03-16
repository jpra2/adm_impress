function [V]=  No_face_element(coord,elem, nsurn1,nsurn2,esurn1,esurn2,nflag,inedge,bedge)
%Retorna os K(n ou t) necessários para a obtenção dos weights. kmap é a
%matriz de permeabilidade; Ex: Kt1->linhaN=Kt1(cellN);


%Verifica quantos são os elementos em torno do nó "ni".%
%V=zeros(2,6,size(coord,1));

for No=1:size(coord,1)
    N_element_No=esurn2(No+1)-esurn2(No);
    
    n_pontos=nsurn2(No+1)-nsurn2(No);
    
    %construção do tensor permeabilidade.%
    
    for k=1:N_element_No
        
        n1=elem(esurn1(esurn2(No)+k),1);
        n2=elem(esurn1(esurn2(No)+k),2);
        n3=elem(esurn1(esurn2(No)+k),3);
        n4=elem(esurn1(esurn2(No)+k),4);
        a=zeros(2,1);
        ii=1;
        for jj=[n1,n2,n3,n4]
            if jj~=No && jj==0
                a(ii,1)=jj;
                ii=ii+1;
            elseif jj~=No && jj~=0
                for g=1:n_pontos
                    h=nsurn1(nsurn2(No)+g);
                    if jj==h
                        a(ii,1)=jj;
                        ii=ii+1;
                    end
                end
            end
        end
        list=nsurn1(nsurn2(No)+1:nsurn2(No+1));
        list2=esurn1(esurn2(No)+1:esurn2(No+1));
        
        for g=1:size(list,1)
            h=list(g);
            if length(list)==length(list2)
                if length(list)==k
                    if a(1,1)==list(g)
                        r=size(list,1)+1;
                    elseif a(2,1)==h
                        s=length(list);
                    end
                else
                    if a(1,1)==h
                        r=g;
                    elseif a(2,1)==h
                        s=g;
                    end
                end
            else
                if a(1,1)==h
                    r=g;
                elseif a(2,1)==h
                    s=g;
                end
            end
        end
        if nflag(No,1)==50000
            
            if r==n_pontos & s==n_pontos
                
                
                m=a(2,1);
                n=a(1,1);
                
            elseif r<s
                m=a(1,1);
                n=a(2,1);
            else
                m=a(2,1);
                n=a(1,1);
            end
            
            
        else
            if r<s
                m=a(1,1);
                n=a(2,1);
            else
                m=a(2,1);
                n=a(1,1);
            end
        end
        
        
        [Tt]=faces_no(bedge, inedge,No,m,n);
        V(:,k,No)= Tt';
        
        
    end
end

end

