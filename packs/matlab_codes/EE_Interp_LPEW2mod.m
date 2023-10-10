function [ E,EE ] =EE_Interp_LPEW2mod(O, P, T, Qo, ni, kmap)

global esurn2

% Retorna os valores dos coeficiente E e Ë

%Prealocação dos vetores%
E=zeros(esurn2(ni+1)-esurn2(ni),2);
EE=zeros(esurn2(ni+1)-esurn2(ni),2);

    %Variáveis%
    %R %matriz de rotação%
    %K %tensor de permeabilidade%
    %distancias %vetor de distancias ortogonais%
    %alfa %Parâmetro que irá variar entre os centros das "células" K,Ksigma e Ktau%
    %beta %Parâmetro que irá variar entre as arestas sigma e tau%	    

    %Para cada elemento, deverá ser calculado 2 valores de E%
    %Para cada elemento, deverá ser calculado 2 valores de Ë%


%Loop que percorre os elementos em torno do nó "ni".%

for k=1:size(O,1)

    %Equação Base%
    %E(k,1:2)=(Eksigma,Ektau)%

    %Preenchimento da segunda coluna do vetor E%
  
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	E(k,2)=dot(K*R*talfa(k,2),(R*talfa(k,2)))/distancias(k,2);  %Esse é o Ektau%

        %Pode-se variar os parâmetros para encontrar as células referentes a Ksigma e Ktau%
        %E(k-1,2)=(dot((K*R*talfa(k,2),(R*talfa(k,2))))/distancias(k,2);%  %Esse é o Eksigma,tau%
        %E(k+1,2)=(dot((K*R*talfa(k,2),(R*talfa(k,2))))/distancias(k,2);%  %Esse é o Ektau,tau%
      
    
    %%%%%%%Fim do Preenchimento da segunda coluna do vetor E %%%%%%%
    
    
    %Preenchimento da primeira coluna do vetor E%
    
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	E(k,1)=dot(K*R*talfa(k,1),(R*talfa(k,1)))/distancias(k,1);  %Esse é o Eksigma%

        %Pode-se variar os parâmetros para encontrar as células referentes a Ksigma e Ktau%
        %E(k-1,1)=(dot((K*R*talfa(k,2),(R*talfa(k,2))))/distancias(k,1);%  %Esse é o Eksigma,sigma%
        %E(k+1,1)=(dot((K*R*talfa(k,2),(R*talfa(k,2))))/distancias(k,1);%  %Esse é o Ektau,sigma%

    %%%%%Fim do Preenchimento da primeira coluna do vetor E %%%%%%%

end

for k=1:size(O,1)

    %Equação Base%
    %Ë(k,1:2)=(Eksigma,Ektau)%

    %Preenchimento da segunda coluna do vetor Ë%
  
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	EE(k,2)=dot(K*R*talfa(k,2),(R*(O(k,:)-Qo)))/distancias(k,2);  %Esse é o Ëktau%

        %Pode-se variar os parâmetros para encontrar as células referentes a Ksigma e Ktau%
        %E(k-1,2)=(dot((K*R*talfa(k,2),(R*(O(k-1,:)-Qo))))/distancias(k,2);%  %Esse é o Ëksigma,tau%
        %E(k+1,2)=(dot((K*R*talfa(k,2),(R*(O(k+1,:)-Qo))))/distancias(k,2);%  %Esse é o Ëktau,tau%
      
    
    %%%%%%%Fim do Preenchimento da segunda coluna do vetor Ë %%%%%%%
    
    
    %Preenchimento da primeira coluna do vetor Ë%
    
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	E(k,1)=dot(K*R*talfa(k,1),(R*(O(k,:)-Qo)))/distancias(k,1);  %Esse é o Eksigma%

        %Pode-se variar os parâmetros para encontrar as células referentes a Ksigma e Ktau%
        %E(k-1,1)=(dot((K*R*talfa(k,2),(R*(O(k-1,:)-Qo))))/distancias(k,1);%  %Esse é o Eksigma,sigma%
        %E(k+1,1)=(dot((K*R*talfa(k,2),(R*(O(k+1,:)-Qo))))/distancias(k,1);%  %Esse é o Ektau,sigma%

    %%%%%Fim do Preenchimento da primeira coluna do vetor Ë %%%%%%%

end



