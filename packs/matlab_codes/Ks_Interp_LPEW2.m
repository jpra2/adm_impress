function [ Kt1, Kt2, Kn1, Kn2 ] = Ks_Interp_LPEW2( O, T, Qo, kmap, ni)

global esurn2 esurn1 elem 
%Retorna os K(n ou t) necessários para a obtenção dos weights. kmap é a
%matriz de permeabilidade; Ex: Kt1->linhaN=Kt1(cellN);

nec=esurn2(ni+1)-esurn2(ni);

Kt1=zeros(nec,2); %As colunas representam i=1 e i=2.
Kt2=zeros(nec,1);
Kn1=zeros(nec,2);
Kn2=zeros(nec,1);
K=zeros(3);
K1=zeros(3);
R=[0 1 0; -1 0 0; 0 0 0];

%Construção do tensor permeabilidade.%

%Cálculo das primeiras constantes, para todas as células que concorrem num%
%nó "ni".                                                                 %
for k=1:nec
 
j=esurn1(esurn2(ni)+k);
    
    for i=1:2
        if (size(T,1)==size(O,1))&&(k==nec)&&(i==2)
            K(1,1)=kmap(elem(j,5),2);
            K(1,2)=kmap(elem(j,5),3);
            K(2,1)=kmap(elem(j,5),4);
            K(2,2)=kmap(elem(j,5),5);
            
            Kn1(k,i)=((R*(T(1,:)-Qo)')'*K*(R*(T(1,:)-Qo)'))/norm(T(1,:)-Qo)^2;
            Kt1(k,i)=((R*(T(1,:)-Qo)')'*K*(T(1,:)-Qo)')/norm(T(1,:)-Qo)^2;
        else
            K(1,1)=kmap(elem(j,5),2);
            K(1,2)=kmap(elem(j,5),3);
            K(2,1)=kmap(elem(j,5),4);
            K(2,2)=kmap(elem(j,5),5);
            
            Kn1(k,i)=((R*(T(k+i-1,:)-Qo)')'*K*(R*(T(k+i-1,:)-Qo)'))/norm(T(k+i-1,:)-Qo)^2;
            Kt1(k,i)=((R*(T(k+i-1,:)-Qo)')'*K*(T(k+i-1,:)-Qo)')/norm(T(k+i-1,:)-Qo)^2;
        end
    end
     
    
    %------------------------- Tensores ----------------------------------%
    
    K1(1,1)= kmap(elem(j,5),2);
    K1(1,2)= kmap(elem(j,5),3);
    K1(2,1)= kmap(elem(j,5),4);
    K1(2,2)= kmap(elem(j,5),5);
    
    if (size(T,1)==size(O,1))&&(k==nec)
        
        %------------ Calculo dos K's internos no elemento ---------------%
    
        Kn2(k)=((R*(T(1,:)-T(k,:))')'*K1*(R*(T(1,:)-T(k,:))'))/norm(T(1,:)-T(k,:))^2;
        Kt2(k)=((R*(T(1,:)-T(k,:))')'*K1*(T(1,:)-T(k,:))')/norm(T(1,:)-T(k,:))^2;
    else
        
        Kn2(k)=(R*(T(k+1,:)-T(k,:))')'*K1*(R*(T(k+1,:)-T(k,:))')/norm(T(k+1,:)-T(k,:))^2;
        Kt2(k)=((R*(T(k+1,:)-T(k,:))')'*K1*(T(k+1,:)-T(k,:))')/norm(T(k+1,:)-T(k,:))^2;
    end

end

end

