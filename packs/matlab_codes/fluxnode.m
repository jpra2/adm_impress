function [Flux ,jelem]= fluxnode(esurn1, esurn2,...
    inedge,no,influx,w,nsurn2,faces_esurn,ielem)

vetor_elem=esurn1(esurn2(no)+1:esurn2(no+1)); % vetor de elementos na vizinhança
% do vértice em questão 'no'
vetor_face= faces_esurn(no,1:(nsurn2(no+1)-nsurn2(no)));%vetor de faces
% na vizinhança do vértice v1

vetor3=vetor_elem==ielem; % indice do vetor que coincide
% com o elemento atual "ielem"

r1=find(vetor3==1); % posição que ocupa o 'ielem' no vetor3
% Alocação do vetor de maneira conveniente tendo em conta como primeiro
% elemento 'ielem' e 'iface', respectivamente.
vetor_elem_ordenado = [vetor_elem(r1:size(vetor3,1),1)' vetor_elem(1:r1-1,1)'];
b=vetor_elem_ordenado(1,2); % elemento adjacente a "ielem"
jelem= vetor_elem_ordenado(1,3); % elemento oposto a "ielem"

%% Alocação das faces de maneira conveniente em sentido antihorario
k=1;
for i=1:length(vetor_face)
    face0 = vetor_face(1,i);
    A1=inedge(face0,3);
    A2=inedge(face0,4);
    
    if (A1==ielem && A2==b) || (A1==b && A2==ielem)
        r=face0; % O "ielem" é unico elemento  compartilhan a
        % "iface" com o elemento "b"
        r2=k; % posição que ocupa o iface no vetor "vetor_face"
    end
    k=k+1;
end

if r==vetor_face(length(vetor_face))
    vetor_face_ordenado=[vetor_face(1, r2) vetor_face(1, 1:r2-1)];
else
    vetor_face_ordenado=[vetor_face(1, r2:size(vetor_face,2)) vetor_face(1, 1:r2-1)];
end

%%
% os dois primeiros elementos serão usados em sentido Antihorario
% e os dois ultimos em sentido Horario em geral.

kk=1;%  inicializa o contador
% loop sobre os elementos do vetor_face_ordenado: este vetor ordena
% todas as faces que concorren nó em questão em sentido
% antihorario.
for face=vetor_face_ordenado % t: devolve o número ou flag da face
    % loop de elementos na vizinhança do vértice em questão.
    if (kk<2 || kk==2)% consideremos fluxos em sentido ANTIHORARIO
        ielem_atual=vetor_elem_ordenado(1,kk);
        
        if (inedge(face,3)==ielem_atual && influx(face)>0)||... % fluxo saindo
                (inedge(face,4)==ielem_atual && influx(face)<0)
            flux_correct= abs(influx(face));
        else
            flux_correct=0;
        end
        
        flux_shapeu(kk)= flux_correct;
        
        kk=kk+1;% contador
        
    else %if (kk==3 || kk==4) % em outro caso fluxos em sentido HORARIO
        
        ielem_atual=vetor_elem_ordenado(kk);
        
        if (inedge(face,3)==ielem_atual && influx(face)<0 ) || ... % o fluxo esta entrado
                (inedge(face,4)==ielem_atual && influx(face)>0 )
            
            flux_correct=abs(influx(face));
            
        else
            
            flux_correct=0;
        end
        
        flux_shapeu(kk)=flux_correct;
        
        kk=kk+1;% contador
    end
    
end
Flux = w*sum(flux_shapeu);
end