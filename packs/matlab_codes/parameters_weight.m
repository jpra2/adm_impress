% Calcula os parametros fisicos e geometricos do artigo  Shuia Miao e
% Jiming Wu, 2021. Os parametros encontra-se na pagina 7 equacao 15.

function [netas,E,EE]=parameters_weight(no,kmap,epsilon)
global nsurn2 nsurn1 esurn2 esurn1 elem coord centelem

% vetor de baricentro na vizinhança do nó "ni"
O=zeros(esurn2(no+1)-esurn2(no),3); 
nec=size(O,1);
K=zeros(3);
% Matriz de rotacao em sentido horario
R=[0 1 0; -1 0 0; 0 0 0];
for k=1:nec
    ielem= esurn1(esurn2(no)+k); % elemento em questão
    
    % calculando vertices para ielem
    vertices=elem(ielem,find(elem(ielem,1:4)~=0));
    
    arroundvertices=nsurn1(nsurn2(no) + 1:nsurn2(no + 1));
    
    if length(arroundvertices)==2
        v=arroundvertices;
    else
        [v]=reordenatevertice(arroundvertices,vertices);
    end
    v1=v(1,1);
    v2=v(2,1);
    % baricentro do elemento em questao
    xilem=centelem(ielem,:);
    % t
    tsigma=coord(v1,:)-coord(no,:);
    ttau=coord(v2,:)-coord(no,:);

    % x sigma tau
    xsigma=coord(no,:)+epsilon*tsigma;
    xtau=coord(no,:)+epsilon*ttau;
    
     vec1=xsigma-coord(no,:);
     vec2=xtau-coord(no,:);

    % normal a face interior da regiao de iteracao
    normal1=(R*-(xtau-xsigma)');
    
    % calculo da area
    sielem=norm(cross(xsigma-coord(no,:),xtau-coord(no,:)))/2;
    
    
    % tensor de permeabilidade
    K(1,1)=kmap(elem(ielem,5),2);
    K(1,2)=kmap(elem(ielem,5),3);
    K(2,1)=kmap(elem(ielem,5),4);
    K(2,2)=kmap(elem(ielem,5),5);
    
    % calculo da distancia ortogonal
    dKsigma=norm(cross(xilem-coord(no,:),tsigma))/norm(tsigma);
    
    dKtau=norm(cross(xilem-coord(no,:),ttau))/norm(ttau);
    
    %% ================================================================
    E(k,1)=dot(K*(R*tsigma'),R*tsigma')/dKsigma;
    E(k,2)=dot(K*(R*ttau'),R*ttau')/dKtau;
    %% ===================================================================
    EE(k,1)=dot(K*(R*tsigma'),R*(xilem-coord(no,:))')/dKsigma;
    
    EE(k,2)=dot(K*(R*ttau'),R*(xilem-coord(no,:))')/dKtau;
    %% ================================================================
    netas(k,1)=dot(K*normal1,(normal1+R*(-vec1')))/(2*sielem);
    
    netas(k,2)=dot(K*normal1,(normal1+R*vec2'))/(2*sielem);
    
end

end
function[v]=reordenatevertice(arroundvertices,vertices)

[lia1,locb1]=ismember(arroundvertices,vertices);
j=1;
v=zeros(2,1);
for i=1:length(locb1)
    if locb1(i)~=0 && i==1 && locb1(length(locb1))~=0
        v(2,1)=arroundvertices(i,1);
        j=1;
    elseif locb1(i)~=0
        v(j,1)=arroundvertices(i,1);
        j=j+1;
    end
end
end