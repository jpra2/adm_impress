function [Flux]=fluxfacecont(inedge,no1,no2,influx,...
    nsurn2,faces_esurn,iface,bflux,ielem)

% Lembrese que quando o fluxo tem valor ngativo, quere dizer que o sentido
% do fluxo é de direita para esquerda. E quando o fluxo tem valor positivo
% o sentido do fluxo é de direita para esquera.
%% nó 1
vetor_face= faces_esurn(no1,1:(nsurn2(no1+1)-nsurn2(no1))); % vetor de faces na vizinhança do vértice v1

%% nó 2
vetor_face2= faces_esurn(no2,1:(nsurn2(no2+1)-nsurn2(no2))); % vetor de faces na vizinhança do vértice v1

%% calcula o fluxo na face em questão para o nó 1

if length(vetor_face)==2 % executa quando o nó esta principalmente na quina
    % do dominio
    % primeira face diferente da face em questão em sentido antihorario
    v1 = vetor_face~=iface+size(inedge,1); % iface +a: número da face atual em vetor_face
    
    r22=find(v1==1); % indice do elemento diferente a face atual
    
    i_bedge=vetor_face(1,r22); % número de face diferente ao número de face atual
    
    i_bedge=i_bedge- size(inedge,1) ;
    
    if bflux(iface)>0 % o fluxo saindo do reservatorio
        
        if bflux(i_bedge)>0
            flux = bflux(i_bedge);
            
        else
            flux = -bflux(i_bedge);
            
        end
        if  (bflux(i_bedge)>0 && flux>0) || (bflux(i_bedge)<0 && flux<0)
            flux=flux;
        else
            flux=0;
        end
    else % o fluxo entrando ao reservatorio
        
        if bflux(i_bedge)>0
            flux = -bflux(i_bedge);
            
        else
            flux = bflux(i_bedge);
        end
        if  (bflux(i_bedge)>0 && flux>0) || (bflux(i_bedge)<0 && flux<0)
            flux=flux;
        else
            flux=0;
        end
    end
    
else % executa quando o nó não esta na quina do dominio
    
    A=vetor_face<=size(inedge,1); % faces interiores
    
    r=find(A==1); % retorna os indices de faces interiores
    
    for j=1:length(r)
        m=vetor_face(r(j)); % escolhe somente faces interiore
        
        a=inedge(m,3);
        b=inedge(m,4);
        
        if (a==ielem || b==ielem) % escolhe, que face interior tem como
            % elemento em comum a ielem.
            face=m;
        end
    end
    
    if bflux(iface)>0 % o fluxo esta saindo do reservatorio
        
        if ielem==inedge(face,4)
            
            
            if influx(face)>0
                fluxo1 = -influx(face);
                
            else
                
                fluxo1 = influx(face);
                
            end
            
        else
            
            if influx(face)>0
                
                fluxo1 = influx(face);
                
            else
                fluxo1 = -influx(face);
                
            end
            
        end
        
        % if ((inedge(face,3)==ielem && influx(face)>0 )|| (inedge(face,4)==ielem && influx(face)<0 )) && ((influx(face)>0 && fluxo1>0) || (influx(face)>0 && fluxo1>0))
        if ((influx(face)>0 && fluxo1>0) || (influx(face)<0 && fluxo1<0))
            flux=fluxo1;
        else
            flux=0;
        end
        
    else % o fluxo entrando ao reservatorio
        
        if ielem==inedge(face,4)
            
            if influx(face)>0
                
                fluxo1 = influx(face);
                
            else
                
                fluxo1 = -influx(face);
                
            end
        else
            if influx(face)>0
                
                fluxo1 = -influx(face);
            else
                fluxo1 = influx(face);
                
            end
            
        end
        
        % if ((inedge(face,3)==ielem && influx(face)>0 )|| (inedge(face,4)==ielem && influx(face)<0 )) && ((influx(face)>0 && fluxo1>0) || (influx(face)>0 && fluxo1>0))
        if  ((influx(face)>0 && fluxo1>0) || (influx(face)<0 && fluxo1<0))
            flux=fluxo1;
        else
            flux=0;
        end
        
    end
    
end

%% Tratamento do fluxo para o nó 2

if length(vetor_face2)==2
    % primeira face diferente da face em questão em sentido antihorario
    v1 = vetor_face2~=iface+size(inedge,1);
    rt= find(v1==1);
    
    i_bedge2=vetor_face2(1,rt);
    
    i_bedge2=i_bedge2 - size(inedge,1);
    
    if bflux(iface)>0
        
        if bflux(i_bedge2)>0
            
            flux1 = bflux(i_bedge2);
            
        else
            flux1 = -bflux(i_bedge2);
            
        end
        if  (bflux(i_bedge2)>0 && flux1>0) || (bflux(i_bedge2)<0 && flux1<0)
            flux1=flux1;
        else
            flux1=0;
        end
    else
        
        if bflux(i_bedge2)>0
            flux1 = -bflux(i_bedge2);
            
        else
            flux1 = bflux(i_bedge2);
            
        end
        if  (bflux(i_bedge2)>0 && flux1>0) || (bflux(i_bedge2)<0 && flux1<0)
            flux1=flux1;
        else
            flux1=0;
        end
    end
    
else
    
    A=vetor_face2<=size(inedge,1);
    r=find(A==1);
    for j=1:length(r)
        m=vetor_face2(r(j));
        
        a=inedge(m,3);
        b=inedge(m,4);
        
        if (a==ielem || b==ielem)
            face2=m;
        end
    end
    
    
    if bflux(iface)>0 % fluxo saindo do reservatorio
        
        if  inedge(face2,3)==ielem
            
            if influx(face2)>0
                
                fluxo2 = influx(face2);
                
            else
                fluxo2 = -influx(face2);
                
            end
        else
            if influx(face2)>0
                
                fluxo2 = -influx(face2);
                
            else
                fluxo2 = influx(face2);
                
            end
            
        end
        
        %if ((inedge(face2,3)==ielem && influx(face2)>0 )|| (inedge(face2,4)==ielem && influx(face2)<0 )) && ((influx(face2)>0 && fluxo2>0) || (influx(face2)>0 && fluxo2>0))
        if ((influx(face2)>0 && fluxo2>0) || (influx(face2)<0 && fluxo2<0))
            flux1=fluxo2;
        else
            flux1=0;
        end
        
    else % fluxo entrando ao reservatorio
        if  inedge(face2,3)==ielem
            
            if influx(face2)>0
                
                fluxo2 = -influx(face2);
                
            else
                fluxo2 = influx(face2);
                
            end
        else
            
            if influx(face2)>0
                
                fluxo2 = influx(face2);
                
            else
                fluxo2 = -influx(face2);
                
            end
            
        end
        
        % if ((inedge(face2,3)==ielem && influx(face2)>0 )|| (inedge(face2,4)==ielem && influx(face2)<0 )) && ((influx(face2)>0 && fluxo2>0) || (influx(face2)>0 && fluxo2>0))
        if ((influx(face2)>0 && fluxo2>0) || (influx(face2)<0 && fluxo2<0))
            
            flux1=fluxo2;
        else
            flux1=0;
        end
    end
    
end

Flux=flux1 + flux;
end