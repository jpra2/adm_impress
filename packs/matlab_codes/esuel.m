function [esuel, esuel_coord,A,bound] = esuel
global elem inedge bedge coord centelem
proj=@(u,v)(dot(u,v)/dot(v,v))*v;
% a ordena��o esta tudo em sentido antihorario 
esuel=zeros(size(elem,1),4);
for i=1:size(elem,1)
    
    % seleciona se triangulo ou quadrilatero
    if elem(i,4)==0
        k=3;
    else
        k=4;
    end
    ll=1;
    for j=1:k % loop de faces
        if j<k
            n1=elem(i,j); % o n� inicio da face en quest�o
            n2=elem(i,j+1); % o n� final da face en quest�o
        else
            n1=elem(i,j); % o n� inicio da face en quest�o
            n2=elem(i,1);% o n� final da face en quest�o
        end
        % calcula o num�ro da face face em quest�o onde est�o contidos o
        % 'n1' e 'n2'.
        a=find((inedge(:,1)==n1 & inedge(:,2)==n2) |(inedge(:,1)==n2 & inedge(:,2)==n1) );
        if a~=0
            if inedge(a,3)~=i % o elemento a esquerda � diferente a elemento atual 'i' fa�a
                esuel(i,j)=inedge(a,3);
            else
                esuel(i,j)=inedge(a,4);
            end
        else
            a=find((bedge(:,1)==n1 & bedge(:,2)==n2) |(bedge(:,1)==n2 & bedge(:,2)==n1) );
            
            esuel(i,j)=bedge(a,3);
            bound(i,ll)=a; % determina o n�mero das faces da fronteira
            ll=ll+1;
            
        end
    end
    
end
% elemento repitidos da matriz 'esuel' somente ocorre quando o elemento
% forma parte do contorno.


% calculando as coordenadas dos elementos fantasma 

for i = 1:size(esuel,1)
    
    s=1;
    for kk=1:size(esuel,2) % loop sobre as colunas da matriz esuel
        
        if esuel(i,kk)~=i & esuel(i,kk)~=0 % evita o indice de um vetor seja zero
            % e excluye as proprierades geometrica do elemento atual 'i'
            
            esuel_coord(kk,:,i)=centelem(esuel(i,kk),:); % Ordena os centroides 
            % dos elementos na vizinhan�a do elemento em quest�o
            
            A(kk,:,i)= centelem(esuel(i,kk),:)-centelem(i,:); % diferen�a de centroides
            
        elseif esuel(i,kk)==i & esuel(i,kk)~=0 % evita o indice de um vetor seja zero
            % e incluye as proprierades geometrica do elemento atual 'i',
            % isto, acontece somente para elemento do contorno.
            for iface =s:size(bound,2)
                if bound(i,iface)~=0
                    
                 % esta pequena rutina calcula o centroide do elemento
                 % fantasma
                coordv1 = coord(bedge(bound(i,iface),1),:);
                coordv2 = coord(bedge(bound(i,iface),2),:);
                xL = centelem(i,:);
                
                p = proj((xL-coordv1),(coordv2-coordv1));
                aux3 = xL-(coordv1+p);
                
                esuel_coord(kk,:,i) = (xL-2*aux3)'; % centroide do elemento
                % fantasma
               
                A(kk,:,i)=(xL-2*aux3)-centelem(i,:); % diferen�a do centroido
                % do elemento 'i' a centroide do elemento fantasma
                s=s+1;
                break
                
                end
            end
        end
        
    end
    
end

end