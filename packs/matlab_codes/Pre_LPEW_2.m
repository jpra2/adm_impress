function [ w,s,wg] = Pre_LPEW_2(kmap,N,gravrate,gravelem)
global coord bcflag bedge nsurn1 nsurn2 gravitational 
% devolve os pesos "w" cujos elementos s�o organizados em um vetor
%Retorna todos os par�metros necess�rios �s express�es dos fluxos.
apw=ones(size(coord,1),1);
apw1=ones(size(coord,1),1);
r=zeros(1,2);
s=zeros(size(bedge,1),1);
global respp;

% kmap(:, 2) = 1;
% kmap(:, 3) = 0;
% kmap(:, 4) = 0;
% kmap(:, 5) = 1;


contg = 1;

id_nodes = 1:size(coord, 1);

for No=1:size(coord,1)
    id_no = id_nodes(No);
    % calcula
    % O--> coordenadas do baricentro na vizinhan�a do n� "No"
    % P--> coordenadas dos v�rtices na vizinhan�a do n� "No"
    % T--> coordenadas dos pontos medios nas vizinhas ao n� "No"
    % Qo-> coordenada do n� em quest�o
    [ O, P, T, Qo ] = OPT_Interp_LPEW(No);
    % calcula os angulos apropiados para calculas os pesos
    
    [ ve2, ve1, theta2, theta1 ] = angulos_Interp_LPEW2( O, P, T, Qo,No );
    % calculas as netas uma relacao de de alturas
    
    [ neta ] = netas_Interp_LPEW( O, P, T, Qo, No );
    % calculas as projecoes normais em torno do n� "No"
    [ Kt1, Kt2, Kn1, Kn2 ] = Ks_Interp_LPEW2( O, T, Qo, kmap, No);
    % calcula os lambdas
    [ lambda,r ] =  Lamdas_Weights_LPEW2( Kt1, Kt2, Kn1, Kn2, theta1,...
        theta2, ve1, ve2, neta, P, O,Qo,No,T,r );
    % gravitational term contribuition
    [g]=gravnode(N,kmap,No,gravelem);
    for k=0:size(O,1)-1
        w(apw(No)+k,1)=lambda(k+1)/sum(lambda); %calculo dos pesos
        respp.coord_nos(contg,:) = Qo;
        respp.coord_faces(contg,:) = O(k+1,:);
        respp.vecw(contg,1) = w(apw(No)+k,1);
        respp.lambdas(contg, 1) = lambda(k+1);
        contg = contg+1;
    end
    
    apw(No+1)=apw(No)+size(O,1);
    
    wg(No,1)=sum(g(:)); %calculo dos pesos
   
    
    %apw1(No+1)=apw1(No)+size(g,1);
    
    % calculando os pesos nos vertices do contorno de Neumann 
    % N ordena faces na vizinhanca de um vertices, comecando pela face do
    % contorno 
    vetor = nsurn1(nsurn2(No) + 1:nsurn2(No + 1));
    comp1 = N(No,1);
    comp2 = N(No,length(vetor));
     MM=bedge(:,1)==No;
     MMM= find(MM == 1);
        if comp1<= size(bedge,1) && comp2 <=size(bedge,1) && 200<bedge(MMM,4)
            % norma das faces
            norm1=norm(coord(bedge(comp1,2),1:2)-coord(bedge(comp1,1),1:2));
            norm2=norm(coord(bedge(comp2,2),1:2)-coord(bedge(comp2,1),1:2));
            a = bcflag(:,1) == bedge(comp1,5);
            s1 = find(a == 1);
            b = bcflag(:,1) == bedge(comp2,5);
            s2 = find(b == 1);
            if strcmp(gravitational,'yes')
                % contribuicoes do termo gravitacional
                m1= gravrate(comp1,1); 
                m2= gravrate(comp2,1);
            else
                m1=0;
                m2=0;
            end
            % da errado quando colocamos o termo gravitacional
            s(No,1) = -(1/sum(lambda))*(r(1,1)*(norm1*bcflag(s1,2) + m1)+...
                                        r(1,2)*(norm2*bcflag(s2,2) + m2));
            %esta rutina funciona quando o vertice da quina da malha
            %computacional pertence ao contorno de Dirichlet
            % s(No,1) = -(1/sum(lambda))*(r(1,1)*bcflag(s1,2) + ...
            %    r(1,2)*bcflag(s2,2));
            %esta rutina funciona quando o vertice da quina da malha
            %computacional pertence ao contorno de Dirichlet ou Neumann
            %s(No,1) = -(1/sum(lambda))*(r(No,1)*m11 + ...
           
            %    r(No,2)*m12);
        end  %End of IF

end
respp.kmap = kmap;
% save('results/resp.mat', '-struct', 'resp');
% disp('saved');
end

