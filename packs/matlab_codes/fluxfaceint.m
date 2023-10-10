function [flux_cros_shapeu]=fluxfaceint(inedge,no1,no2,influx,...
    bedge,w,nsurn2,faces_esurn,iface,bflux,ielem1,jelem1)
% Define o elemento atual


% Lembrese que quando o fluxo tem valor ngativo, quere dizer que o sentido
% do fluxo é de direita para esquerda. E quando o fluxo tem valor positivo
% o sentido do fluxo é de direita para esquera.
vetor_face= faces_esurn(no1,1:(nsurn2(no1+1)-nsurn2(no1))); % vetor de faces na vizinhança do vértice v1

vetor31=vetor_face'==iface; % indice do vetor que coincide
% com a face "iface"

r2=find(vetor31==1); % posição que ocupa a iface em questão no vetor31

% alocação do vetor de maneira conveniente tendo em conta como primeiro
% elemento "ielem" e "iface"

vetor_face_ordenado=[vetor_face(1, r2:size(vetor31,1)) vetor_face(1, 1:r2-1)];
%%
% do vértice v1
vetor_face2= faces_esurn(no2,1:(nsurn2(no2+1)-nsurn2(no2))); % vetor de faces na vizinhança do vértice v1

% com o elemento "ielem"
vetor312=vetor_face2'==iface; % indice do vetor que coincide
% com a face "iface"
r22=find(vetor312==1); % posição que ocupa a iface em questão no vetor31

% alocação do vetor de maneira conveniente tendo em conta como primeiro
% elemento "ielem" e "iface"

vetor_face_ordenado2=[vetor_face2(1, r22:size(vetor312,1)) vetor_face2(1, 1:r22-1)];

%% calcula o fluxo na face em questão

%% no 1

% primeira face diferente da face em questão em sentido antihorario
n_face = vetor_face_ordenado(1,2) ;

% ultima face diferente da face em questão em sentido antihorario
m_face = vetor_face_ordenado(1,size(vetor_face_ordenado,2));

%% Tratamento para o nó 2
n_face2 = vetor_face_ordenado2(1,2) ;

% ultima face diferente da face em questão em sentido antihorario
m_face2 = vetor_face_ordenado2(1,size(vetor_face_ordenado2,2));

tfaces=[n_face m_face n_face2 m_face2];

if size(inedge,1)<n_face
    a=bedge(n_face-size(inedge,1) ,3)==ielem1;
    c=max(a);
    if c~=0
        n_face=n_face;
        m_face=m_face;
    else
        n_face=m_face;
        m_face=tfaces(1);
    end
else
    a=inedge(n_face,3:4)==ielem1;
    c=max(a);
    if c~=0
        n_face=n_face;
        m_face=m_face;
    else
        n_face=m_face;
        m_face=tfaces(1);
    end
end
if  size(inedge,1)<m_face2
    a1=bedge(m_face2-size(inedge,1),3)==jelem1;
    c1=max(a1);
    if c1~=0
        n_face2=n_face2;
        m_face2=m_face2;
    else
        n_face2=m_face2;
        m_face2=tfaces(3);
    end
else
    a1=inedge(m_face2,3:4)==jelem1;
    c1=max(a1);
    if c1~=0
        n_face2=n_face2;
        m_face2=m_face2;
    else
        n_face2=m_face2;
        m_face2=tfaces(3);
    end
end

k=1;
for face=[n_face n_face2]
    
    if size(inedge,1)<face
        i_bedge1=face - size(inedge,1); % face do contorno
        
        if (bflux(i_bedge1)>0 && ielem1==bedge(i_bedge1,3))||(bflux(i_bedge1)<0 && ielem1==bedge(i_bedge1,3))  % fluxo  entrando ao elemento "jelem"
            
            
            flux_face =abs(bflux(i_bedge1));
            
        else
            flux_face = 0;
            
        end
        
        flux_shapeu(k)=flux_face;
        k=k+1;
       
        
    else
        
        
        if (influx(face)<0 && inedge(face,4)==ielem1)||(influx(face)>0 && inedge(face,3)==ielem1) % fluxo entrando a "jelem"
            flux=abs(influx(face));
            
        else % fluxo saindo de "ielem: elemento atual"
            flux=0;
            
        end
        flux_shapeu(k)=flux;
        k=k+1;
        
    end
    
end

for face=[m_face m_face2]
    
    if size(inedge,1)<face
        i_bedge2=face - size(inedge,1);
        if (bflux(i_bedge2)>0 && jelem1==bedge(i_bedge2,3))||(bflux(i_bedge2)<0 && jelem1==bedge(i_bedge2,3))  % fluxo  entrando ao elemento "jelem"
            
            
            flux_face =abs(bflux(i_bedge2));
            
        else
            flux_face = 0;
            
        end
        
        flux_shapeu(k)=flux_face;
        k=k+1;
    else
        if (influx(face)>0 && inedge(face,4)==jelem1)||(influx(face)<0 && inedge(face,3)==jelem1)
            flux=abs(influx(face));
        else
            flux=0;
        end
        flux_shapeu(k)=flux;
        k=k+1;
        
    end
end
if (inedge(iface,3)==ielem1 && influx(iface)>0)||(inedge(iface,4)==ielem1 && influx(iface)<0)
    flux_face1=abs(influx(iface));
   
else
    flux_face1=0;
   
end
flux_cros_shapeu=(1-4*w)*flux_face1 + w*sum(flux_shapeu);
end