function [QI,QJ] = Q_elem(coord,esurn1,esurn2,c,xD,xM,inedge)

QI = zeros(size(inedge,1),1);
QJ = zeros(size(inedge,1),1);

% the angle between two vectors is the "dot" product divided by the product
% of modules
angle = @(u,v)acos(dot(u,v)/(norm(u)*norm(v)));

for i = 1:size(inedge,1)
    
    angI = 1e10;
    angJ = 1e10;
    
    v1 = inedge(i,1);
    v2 = inedge(i,2);
    
    lef = inedge(i,3); % elemento a esquerda
    rel = inedge(i,4); % elemento a direita
    
    if norm(coord(v2,:)-xD(i,:)) > (norm(coord(v2,:)-coord(v1,:))/2)
        P = v2;
    else
        P = v1;
    end
    
    rDM = xM(i,:)-xD(i,:); % reta DM
    rrr=norm(rDM);
    
    % retorna uma lista de elementos na vizinhança do vértice P
    
    list = esurn1((esurn2(P)+1):esurn2(P+1));
    
    while list(1)~=lef 
        aux = list(1);
        for j = 2:size(list,1)
            list(j-1) = list(j);
        end
        list(size(list,1)) = aux;
    end
    
    for j = 2:size(list,1)
        elem = list(j);
        r = c(elem,:)-c(lef,:);
        % angle = @(u,v)acos(dot(u,v)/(norm(u)*norm(v)));
        %----------------------------------------------------------
        if norm(rDM) < 1e-10
            
            ang1=angle(r,(coord(v2,:)-coord(v1,:)));
            
        else
            
            ang1 = angle(r,rDM);
        end
        %----------------------------------------------------------
         %ang1 = angle(r,rDM);
        % calcula o menor angulo
        if ang1 < angI
            angI = ang1;
            QI(i) = elem;
        end
    end
    
    
    while list(1)~=rel
        aux = list(1);
        for j = 2:size(list,1)
            list(j-1) = list(j);
        end
        list(size(list,1)) = aux;
    end
    
    for j = 2:size(list,1)
        elem = list(j);
        r = (c(elem,:)-c(rel,:))';
        %---------------------------------------------------
        if norm(rDM) < 1e-10
            
            ang2=angle(r,(coord(v2,:)-coord(v1,:)));
            
        else
            ang2 = angle(r,rDM);
        end
        %---------------------------------------------------
         %ang2 = angle(r,rDM);
        % calcula o menor angulo
        
        if ang2 < angJ
            angJ = ang2;
            QJ(i) = elem;
        end
    end
end
end