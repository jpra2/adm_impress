% objetivo: Montagem da matriz global M e I
function [ M, I ] = globalmatrixmpfad( w,s, Kde, Ded, Kn, Kt, nflagno, ...
    Hesq,fonte,gravresult,gravrate,gravno,gravelem,grav_elem_escalar,wg)

global coord elem esurn1 esurn2  bedge inedge  centelem bcflag gravitational...
    strategy elemarea

global respp;

cont_dir = 1;
cont_neumann = 1;
cont_neumann_edges = 1;
cont_neumann_weight = 1;
dirichlet_nodes = zeros(1,3);
neumann_nodes = zeros(1,3);
neumann_edges = zeros(1,3);
neumann_weight = zeros(1,3);

%-----------------------inicio da rutina ----------------------------------%
%Constroi a matriz global.

M=sparse(size(elem,1),size(elem,1)); %Prealocacao de M.
I=sparse(size(elem,1),1);
% contribuicao do termo de fonte
% I=I+fonte;

m=0;
cont_edge = 1;
cont_face = 1;

for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:); %fase. b1b2
    v1=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,1),:); % b1ob
    v2=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,2),:); % b2ob
    edge_centroid = (coord(bedge(ifacont,2),:) + coord(bedge(ifacont,1),:))./2;
    edges_centroids(cont_edge, :) = edge_centroid(1:2);
    normcont=norm(v0);
    
    notest = 100;
    n1 = bedge(ifacont,1);
    n2 = bedge(ifacont,2);
    
    if n1 == notest | n2 == notest
        disp('');
    end
    
    
    % Tratamento do n� nos v�rtices 2 e 4%
    A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
    ked_data(cont_edge,:) = A;
    kn_data(cont_edge,:) = Kn(ifacont);
    kt_data(cont_edge,:) = Kt(ifacont);
    ded_data(cont_edge,:) = 0;
    
    cont_edge = cont_edge + 1;
    
    if bedge(ifacont,5)<200
        % Contorno de Dirichlet
        c1=nflagno(bedge(ifacont,1),2);
        c2=nflagno(bedge(ifacont,2),2);
        
        %Preenchimento do termo gravitacional
        
        if strcmp(gravitational,'yes')
            if strcmp(strategy,'starnoni')||strcmp(strategy,'inhouse3')
                m=gravrate(ifacont);
            elseif strcmp(strategy,'inhouse1')
                g1=gravno(bedge(ifacont,1),1); % gravidade no vertice 1
                g2=gravno(bedge(ifacont,2),1); % gravidade no vertice 2
               % m=-(A*(dot(v2,-v0)*g1+dot(v1,v0)*g2-(norm(v0)^2*grav_elem_escalar(lef)))-(g2-g1)*Kt(ifacont));
                m=gravrate(ifacont)+g2+g1;
            end
        else
            m=0;
        end
        % ambos os vertices pertenecem ao contorno de Dirichlet
        if nflagno(bedge(ifacont,2),1)<200 && nflagno(bedge(ifacont,1),1)<200
            if bedge(ifacont,2) == 100 | bedge(ifacont,1)==100
                disp('');
            end
            %montagem da matriz global 
            M(lef,lef)=M(lef,lef)-A*(norm(v0)^2);
            % termo de fonte
            I(lef)=I(lef)-A*(dot(v2,-v0)*c1+dot(v1,v0)*c2)+(c2-c1)*Kt(ifacont)+m;
            
            data_to_insert(1, 1:2) = coord(bedge(ifacont,1),1:2);
            data_to_insert(1, 3) = c1;
            [dirichlet_nodes, cont_dir] = insert_data_list(dirichlet_nodes, cont_dir, data_to_insert);
            data_to_insert(1, 1:2) = coord(bedge(ifacont,2),1:2);
            data_to_insert(1, 3) = c2;
            [dirichlet_nodes, cont_dir] = insert_data_list(dirichlet_nodes, cont_dir, data_to_insert);
            
        else
            % quando um dos vertices da quina da malha computacional
            % pertence ao contorno de Neumann
            if nflagno(bedge(ifacont,1),1)>200
                %montagem da matriz global
                M(lef,lef)=M(lef,lef)-A*(norm(v0)^2)+Kt(ifacont)+A*dot(v2,-v0);
                % termo de fonte
                I(lef)=I(lef)-A*(dot(v1,v0)*c2)+(c2)*Kt(ifacont)+m;
                
                data_to_insert(1, 1:2) = coord(bedge(ifacont,2),1:2);
                data_to_insert(1, 3) = c2;
                [dirichlet_nodes, cont_dir] = insert_data_list(dirichlet_nodes, cont_dir, data_to_insert);
                
                
            elseif nflagno(bedge(ifacont,2),1)>200
                %montagem da matriz global
                M(lef,lef)=M(lef,lef)-A*(norm(v0)^2)-Kt(ifacont)+A*dot(v1,v0);
                % termo de fonte
                I(lef)=I(lef)-A*(dot(v2,-v0)*c1)+(-c1)*Kt(ifacont)+m;
                
                data_to_insert(1, 1:2) = coord(bedge(ifacont,1),1:2);
                data_to_insert(1, 3) = c1;
                [dirichlet_nodes, cont_dir] = insert_data_list(dirichlet_nodes, cont_dir, data_to_insert);
                disp(data_to_insert(3));
                disp('');
            end
        end
    else
        % Contorno de Neumann
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef) -normcont*bcflag(r,2);
        
        data_to_insert(1, 1:2) = edge_centroid(1:2);
        data_to_insert(1, 3) = bcflag(r,2);
        [neumann_edges, cont_neumann_edges] = insert_data_list(neumann_edges, cont_neumann_edges, data_to_insert);
    end 
end

% contribui��o nas faces internas
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    vol_test = 1;
    if lef == vol_test | rel == vol_test
        disp('');
    end
    
    edge_centroid = (coord(inedge(iface,2),:) + coord(inedge(iface,1),:))./2;
    edges_centroids(cont_edge, :) = edge_centroid(1:2);
    
    ked_data(cont_edge,:) = Kde(iface);
    kn_data(cont_edge,:) = 0;
    kt_data(cont_edge,:) = 0;
    ded_data(cont_edge,:) = Ded(iface);
    
    cont_edge = cont_edge + 1;
    
    %Contabiliza as contribui��es do fluxo numa aresta para os elementos %
    %a direita e a esquerda dela.                                        %
    
    M(lef, lef)=M(lef, lef)- Kde(iface);
    M(lef, rel)=M(lef, rel)+ Kde(iface);
    M(rel, rel)=M(rel, rel)- Kde(iface);
    M(rel, lef)=M(rel, lef)+ Kde(iface);
    
    %Se os n�s das arestas estiverem em fronteiras de Dirichlet, suas
    %contribui��es ser�o contabilizadas logo abaixo.
    notest = 100;
    n1 = inedge(iface,1);
    n2 = inedge(iface,2);
    if n1 == notest | n2 == notest
        disp('');
    end
    
    
    if nflagno(inedge(iface,1),1)<200
        I(lef)=I(lef)-Kde(iface)*Ded(iface)*nflagno(inedge(iface,1),2);
        I(rel)=I(rel)+Kde(iface)*Ded(iface)*nflagno(inedge(iface,1),2);
    end
    if nflagno(inedge(iface,2),1)<200
        I(lef)=I(lef)+Kde(iface)*Ded(iface)*nflagno(inedge(iface,2),2);
        I(rel)=I(rel)-Kde(iface)*Ded(iface)*nflagno(inedge(iface,2),2);
    end
    % quando o n� pertece ao contorno de Neumann
    if nflagno(inedge(iface,1),1)==201
        
        I(inedge(iface,3))=I(inedge(iface,3))-Kde(iface)*Ded(iface)*s(inedge(iface,1)); %ok
        
        I(inedge(iface,4))=I(inedge(iface,4))+Kde(iface)*Ded(iface)*s(inedge(iface,1)); %ok
        
        data_to_insert(1, 1:2) = coord(bedge(ifacont,1),1:2);
        data_to_insert(1, 3) = s(inedge(iface,1));
        [neumann_weight, cont_neumann_weight] = insert_data_list(dirichlet_nodes, cont_neumann_weight, data_to_insert);
            
    end
    if nflagno(inedge(iface,2),1)==201
        
        I(inedge(iface,3))=I(inedge(iface,3))+Kde(iface)*Ded(iface)*s(inedge(iface,2)); %ok
        
        I(inedge(iface,4))=I(inedge(iface,4))-Kde(iface)*Ded(iface)*s(inedge(iface,2)); %ok
        
        data_to_insert(1, 1:2) = coord(bedge(ifacont,2),1:2);
        data_to_insert(1, 3) = s(inedge(iface,2));
        [neumann_weight, cont_neumann_weight] = insert_data_list(dirichlet_nodes, cont_neumann_weight, data_to_insert);
        
        
    end
    
    %Contabiliza��o das contribui��es dos n�s que n�o est�o na
    %fronteiras de Dirichlet.
    % first node
    if nflagno(inedge(iface,1),1)>200
        for j=1:(esurn2(inedge(iface,1)+1)-esurn2(inedge(iface,1)))
            
            post_cont=esurn2(inedge(iface,1))+j;
            
            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
            
            M(rel, esurn1(post_cont))=M(rel,esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
           
        end
        %I(lef)=I(lef)+  wg(inedge(iface,1));
        %I(rel)=I(rel)-   wg(inedge(iface,1));
    end
    % second node
    if nflagno(inedge(iface,2),1)>200
        for j=1:(esurn2(inedge(iface,2)+1)-esurn2(inedge(iface,2)))
            
            post_cont=esurn2(inedge(iface,2))+j;
            
            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
            
            M(rel, esurn1(post_cont))=M(rel,esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
            
        end
        %I(lef)=I(lef)-  wg(inedge(iface,2));
        %I(rel)=I(rel)+  wg(inedge(iface,2));
    end
    % termo gravitacional
    if strcmp(gravitational,'yes')
        if strcmp(strategy,'starnoni') ||strcmp(strategy,'inhouse3')
            m=gravrate(size(bedge,1)+iface,1);
        elseif strcmp(strategy,'inhouse')
            no1=inedge(iface,1);
            no2=inedge(iface,2);
            %g1=gravno(no1,1);
            %g2=gravno(no2,1);
            g1=0;
            nec1=esurn2(no1+1)- esurn2(no1);
            nec2=esurn2(no2+1)- esurn2(no2);
            for j=1:nec1
                    element1=esurn1(esurn2(no1)+j);
                    g1=g1+w(esurn2(no1)+j)*grav_elem_escalar(element1);
            end
            g2=0;
            for j=1:nec2
                    element2=esurn1(esurn2(no2)+j);
                    g2=g2+w(esurn2(no2)+j)*grav_elem_escalar(element2);
            end
           % m= -Kde(iface)*(grav_elem_escalar(rel)-grav_elem_escalar(lef)-Ded(iface)*(g2-g1));
            m= gravrate(size(bedge,1)+iface,1)+g2+g1;
        else
            m=0;
        end
        
        I(lef)=I(lef)+m ;
        I(rel)=I(rel)-m ;
    end
end
[rowM, colM, dataM] = find(M);
respp.rowM = rowM;
respp.colM = colM;
respp.dataM = dataM;
respp.source = full(I);

respp.dirichlet_nodes = dirichlet_nodes;
respp.neumann_nodes = neumann_nodes;
respp.neumann_edges = neumann_edges;
respp.neumann_weight = neumann_weight;

respp.centelem = centelem;
respp.coord = coord;
respp.fonte = fonte;
respp.edges_centroids = edges_centroids;

respp.ked_data = ked_data;
respp.kn_data = kn_data;
respp.kt_data = kt_data;
respp.ded_data = ded_data;


I = I+fonte;
end
