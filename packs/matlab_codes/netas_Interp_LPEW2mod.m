function [ netas ] =netas_Interp_LPEW2mod(O, P, T, Qo, ni, kmap)

global esurn2
% Retorna os netas.

%Prealocação do vetor.%
netas=zeros(esurn2(ni+1)-esurn2(ni),2);


%Variáveis de interesse%
    %nkv %normal relativa a xsigma-xtau%
    %nksigma %normal relativa a xv-xtau%
    %nktau %normal relativa a xv-xsigma%
    %K  %Tensor de permeabilidade%
    %R  %Matriz de Rotação%
    %Skv %area%
    %beta %Parâmetro que irá variar entre as arestas%
    %sigma %Primeira Aresta%
    %tau %Segunda Aresta%	    

    %Para cada elemento, deverá ser calculado 2 valores de neta%
    %Equação Base%
    %netas(k,1)=(dot((K*nkv),(nkv*norm(beta)*nksigma)))/2*Skv%
    %netas(k,1:2)=(netasigma,netatau)%

%Loop que percorre os elementos em torno do nó "ni".%

for k=1:size(netas,1)


    %Preenchimento da segunda coluna do vetor "netas"%
    
    if (k==size(netas,1))&&(size(P,1)==size(O,1))
        tau=P(1,:)-Qo;
        vetor=T(2,:)-T(1,:);
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	nkv=(R*vetor)/(norm(R*vetor));
	netas(k,2)=(dot((K*nkv),(nkv*norm(tau)*nktau)))/2*Skv; %Este é o neta,ktau%

    else
        tau=P(k+1,:)-Qo;
	vetor=T(k+1,:)-T(k,:);
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	nkv=(R*vetor)/(norm(R*vetor));
	netas(k,2)=(dot((K*nkv),(nkv*norm(tau)*nktau)))/2*Skv; %Este é o neta,ktau%

    end
    
    %%%%%%%Fim do Preenchimento da segunda coluna do vetor "netas".%%%%%%%
    
    
    %Preenchimento da primeira coluna do vetor "netas".%
    
        sigma=P(1,:)-Qo;
        vetor=T(2,:)-T(1,:);
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	nkv=(R*vetor)/(norm(R*vetor));
	netas(k,1)=(dot((K*nkv),(nkv*norm(sigma)*nksigma)))/2*Skv; %Este é o neta,ksigma%
        
    %%%%%Fim do Preenchimento da primeira coluna do vetor "netas".%%%%%%%

end

