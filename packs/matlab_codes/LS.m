function [ w,s] = LS(kmap)
global coord bcflag bedge nsurn1 nsurn2 esurn1 esurn2 inedge elem centelem normals
% devolve os pesos "w" cujos elementos são organizados em um vetor
%Retorna todos os parâmetros necessários às expressões dos fluxos.
multvec = @(x,k) [x(1)*k(1) + x(2)*k(3),x(1)*k(2) + x(2)*k(4)];
w = 0;
s = 0;
ang =  pi/2;
R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
[bc_matrix] = bcGLS();
w = [];
s = [];
for node = 1:size(coord,1)
    around_elem = esurnG(node);
    around_node = nsurnG(node);
    center_node_type = bc_matrix(node,1);
    start_node_type =  bc_matrix(around_node(1),1);
    end_node_type =  bc_matrix(around_node(end),1);
    start_flux =  bc_matrix(node,3);
    end_flux =  bc_matrix(node,4);
    start_edge =   bc_matrix(node,end-1);
    end_edge =   bc_matrix(node,end);
    %around_k = zeros(size(around_elem,1),4);
    kmap_loc = elem(around_elem,end);
    around_k = kmap(kmap_loc,2:end);
    node_center = coord(node,1:2);
    
    around_elemcenter = centelem(around_elem,1:2);
    around_nodecenter = coord(around_node,1:2);
    
    % calculating tal and tal unitary normals
    tal_int = repelem(node_center ,size(around_nodecenter ,1),1);
    tal = around_nodecenter  -  tal_int;
    length_tal =  sum(tal.^2,2).^.5;
    tal_norm = [R(1,1)*tal(:,1) + R(1,2)*tal(:,2), R(2,1)*tal(:,1) + R(2,2)*tal(:,2)];
    tal_versor = [tal_norm(:,1)./length_tal tal_norm(:,2)./length_tal];
    % calculating xnode - xn
    tal_int = repelem(node_center ,size(around_elemcenter ,1),1);
    xdiff = around_elemcenter - tal_int;
    % calculating tal_Versor times K
    B1 = zeros(size(around_elemcenter,1), size(around_elemcenter,1)*2 + 1);
    jj = 1;
    for ii = 1:size(xdiff,1)
        vec = xdiff(ii,:);
        B1(ii,jj) = vec(1);
        B1(ii,jj+1) = vec(2);
        jj = jj + 2;
    end
    
    B1(:,end) = 1;
    B2 = zeros(size(around_nodecenter,1), size(around_elemcenter,1)*2 + 1);
    
    
    if center_node_type == 0
        jj = 1;
        for ii = 2:size(tal,1)
            vec = tal(ii,:);
            B2(ii,jj) = -vec(1);
            B2(ii,jj+1) = -vec(2);
            B2(ii,jj+2) = vec(1);
            B2(ii,jj+3) = vec(2);
            jj = jj + 2;
        end
        B2(1,1:2) = tal(1,:);
        B2(1,end-2:end-1) = -tal(1,:);
    elseif center_node_type == 1
        1;
        %disp("B2 Dirichlet node");
    elseif center_node_type == 2
        %disp("B2 Neumman Node")
        jj = 1;
        for ii = 2:(size(tal,1)-1)
            vec = tal(ii,:);
            B2(ii,jj) = -vec(1);
            B2(ii,jj+1) = -vec(2);
            B2(ii,jj+2) = vec(1);
            B2(ii,jj+3) = vec(2);
            jj = jj + 2;
        end
        B2(1,1:2) = tal(1,:);
        B2(end,end-2:end-1) = -tal(end,:);
    end
    B2(1,:) = [];
    B2(end,:) = [];
    %%%%
    B3 = zeros(size(around_k ,1), size(around_elemcenter,1)*2 + 1);
    if center_node_type == 0
        B3 = zeros(size(around_k ,1), size(around_elemcenter,1)*2 + 1);
        jj = 1;
        for ii = 2:size(around_k,1)
            vec = tal_versor(ii,:);
            K1 = around_k(ii-1,:);
            K2 = around_k(ii,:);
            vec1 = multvec(vec,K1);
            vec2 = multvec(vec,K2);
            B3(ii,jj) =   -vec1(1);
            B3(ii,jj+1) = -vec1(2);
            B3(ii,jj+2) =  vec2(1);
            B3(ii,jj+3) =  vec2(2);
            jj = jj + 2;
        end
        vec = tal_versor(1,:);
        K1 = around_k(1,:);
        K2 = around_k(end,:);
        vec1 = multvec(vec,K1);
        vec2 = multvec(vec,K2);
        B3(1,1:2) =   -vec1;
        B3(1,end-2:end-1) =  vec2;
    elseif  center_node_type == 1
        1;
        %disp("B3 Dirichlet node")
    elseif  center_node_type == 2
        B3 = zeros(size(around_k ,1) + 1, size(around_elemcenter,1)*2 + 1);
        jj = 1;
        start_normal = normals(start_edge,1:2);
        end_normal = normals(end_edge,1:2);
        start_talnormal = tal_versor(1,:);
        end_talnormal = tal_versor(2,:);
        
        
        for ii = 2:size(around_k,1)-1
            vec = tal_versor(ii,:);
            K1 = around_k(ii-1,:);
            K2 = around_k(ii,:);
            vec1 = multvec(vec,K1);
            vec2 = multvec(vec,K2);
            B3(ii,jj) =   -vec1(1);
            B3(ii,jj+1) = -vec1(2);
            B3(ii,jj+2) =  vec2(1);
            B3(ii,jj+3) =  vec2(2);
            jj = jj + 2;
        end
        %first
        vec = tal_versor(1,:);
        K1 = around_k(1,:);
        %K2 = around_k(end,:);
        vec1 = multvec(vec,K1);
        %vec2 = multvec(vec,K2);
        B3(1,1:2) =   vec1;
        %B3(1,end-2:end-1) =  vec2;
        %last
        vec = tal_versor(end,:);
        K1 = around_k(end,:);
        vec1 = multvec(vec,K1);
        %vec2 = multvec(vec,K2);
        B3(end,end-2:end-1) =   -vec1;
        
        
    end
    M = [B1; B2; B3];
    f = zeros(size(B3,1),1);
    f(1) = start_flux;
    f(end) = end_flux;
    F = [zeros(size(B1,1),1); zeros(size(B2,1),1); f] ;
    
    [w1,s1] =  solvMN(M, size(around_elem,1),F);
    w = [w w1];
    if center_node_type == 2
        s = [s; s1];
    end
    %     nline = size(M,1);
    
end

%N = [eye(nline); zeros(size(

end

function [w, s] = solvMN(M, naround_elem,F)
[l,c] = size(M);
A = eye(naround_elem);
B = zeros(l-naround_elem, naround_elem);
N = [A;B];
%G = (M'*N)\(M'*M);
%G = inv(M'*M)*(M'*N);
G = (M'*M)\(M'*N);
H = (M'*M)\(M'*F);
w = G(end,:);
s = H(end,:);
end
