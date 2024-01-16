function [ E,EE ] =EE_Interp_LPEW2mod(O, P, T, Qo, ni, kmap)

global esurn2

% Retorna os valores dos coeficiente E e �

%Prealoca��o dos vetores%
E=zeros(esurn2(ni+1)-esurn2(ni),2);
EE=zeros(esurn2(ni+1)-esurn2(ni),2);

    %Vari�veis%
    %R %matriz de rota��o%
    %K %tensor de permeabilidade%
    %distancias %vetor de distancias ortogonais%
    %alfa %Par�metro que ir� variar entre os centros das "c�lulas" K,Ksigma e Ktau%
    %beta %Par�metro que ir� variar entre as arestas sigma e tau%	    

    %Para cada elemento, dever� ser calculado 2 valores de E%
    %Para cada elemento, dever� ser calculado 2 valores de �%


%Loop que percorre os elementos em torno do n� "ni".%

for k=1:size(O,1)

    %Equa��o Base%
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
	E(k,2)=dot(K*R*talfa(k,2),(R*talfa(k,2)))/distancias(k,2);  %Esse � o Ektau%

        %Pode-se variar os par�metros para encontrar as c�lulas referentes a Ksigma e Ktau%
        %E(k-1,2)=(dot((K*R*talfa(k,2),(R*talfa(k,2))))/distancias(k,2);%  %Esse � o Eksigma,tau%
        %E(k+1,2)=(dot((K*R*talfa(k,2),(R*talfa(k,2))))/distancias(k,2);%  %Esse � o Ektau,tau%
      
    
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
	E(k,1)=dot(K*R*talfa(k,1),(R*talfa(k,1)))/distancias(k,1);  %Esse � o Eksigma%

        %Pode-se variar os par�metros para encontrar as c�lulas referentes a Ksigma e Ktau%
        %E(k-1,1)=(dot((K*R*talfa(k,2),(R*talfa(k,2))))/distancias(k,1);%  %Esse � o Eksigma,sigma%
        %E(k+1,1)=(dot((K*R*talfa(k,2),(R*talfa(k,2))))/distancias(k,1);%  %Esse � o Ektau,sigma%

    %%%%%Fim do Preenchimento da primeira coluna do vetor E %%%%%%%

end

for k=1:size(O,1)

    %Equa��o Base%
    %�(k,1:2)=(Eksigma,Ektau)%

    %Preenchimento da segunda coluna do vetor �%
  
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	EE(k,2)=dot(K*R*talfa(k,2),(R*(O(k,:)-Qo)))/distancias(k,2);  %Esse � o �ktau%

        %Pode-se variar os par�metros para encontrar as c�lulas referentes a Ksigma e Ktau%
        %E(k-1,2)=(dot((K*R*talfa(k,2),(R*(O(k-1,:)-Qo))))/distancias(k,2);%  %Esse � o �ksigma,tau%
        %E(k+1,2)=(dot((K*R*talfa(k,2),(R*(O(k+1,:)-Qo))))/distancias(k,2);%  %Esse � o �ktau,tau%
      
    
    %%%%%%%Fim do Preenchimento da segunda coluna do vetor � %%%%%%%
    
    
    %Preenchimento da primeira coluna do vetor �%
    
	K(1,1)=kmap(elem(k,5),2);
        K(1,2)=kmap(elem(k,5),3);
        K(2,1)=kmap(elem(k,5),4);
        K(2,2)=kmap(elem(k,5),5);
        R(1,1)=0;
        R(1,2)=1;
        R(2,1)=-1;
        R(2,2)=0;
	E(k,1)=dot(K*R*talfa(k,1),(R*(O(k,:)-Qo)))/distancias(k,1);  %Esse � o Eksigma%

        %Pode-se variar os par�metros para encontrar as c�lulas referentes a Ksigma e Ktau%
        %E(k-1,1)=(dot((K*R*talfa(k,2),(R*(O(k-1,:)-Qo))))/distancias(k,1);%  %Esse � o Eksigma,sigma%
        %E(k+1,1)=(dot((K*R*talfa(k,2),(R*(O(k+1,:)-Qo))))/distancias(k,1);%  %Esse � o Ektau,sigma%

    %%%%%Fim do Preenchimento da primeira coluna do vetor � %%%%%%%

end



