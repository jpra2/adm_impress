%funçao que calcula: permeabilidades relativas - Krw e Kro
%                    mobilidades - lambda_w e lambda_o
%                    fluxo fracional - f

function [mobility] = mobilityface(S_old,nw,no,auxflag,S_cont,simu)

global visc satlimit elem inedge esurn1 esurn2 bedge elemarea
if strcmp(simu,'monofasico')
   mobility=ones(size(inedge,1)+size(bedge,1),1); 
else
    mobility=zeros(size(inedge,1)+size(bedge,1),1);
    for i = 1:size(inedge,1) % loop de fases internas para calcular as mobilidade e fluxo fracional
        v1=inedge(i,1);
        v2=inedge(i,2);
        
        n1=zeros(size(elem,1),1);
        A=zeros(size(elem,1),1);
        n2=zeros(size(elem,1),1);
        A1=zeros(size(elem,1),1);
        for j=esurn2(v1)+1:esurn2(v1+1)
            tt=esurn1(j);
            
            Krw1 = ((S_old(tt) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
            
            Kro1 = ((1 - S_old(tt) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
            
            L11 = Krw1/visc(1) +  Kro1/visc(2); % mobilidade total no elemento esquerdo
            
            A(tt)=elemarea(tt);
            
            n1(tt)=A(tt)*L11;
            
        end
        
        L1=sum(n1)/sum(A);
        
        for jj=esurn2(v2)+1:esurn2(v2+1)
            tt=esurn1(jj);
            
            
            Krw1 = ((S_old(tt) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
            
            Kro1 = ((1 - S_old(tt) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
            
            L22 = Krw1/visc(1) +  Kro1/visc(2); % mobilidade total no elemento esquerdo
            
            A1(tt)=elemarea(tt);
            
            n2(tt)=A1(tt)*L22;
            
        end
        
        L2=sum(n2)/sum(A1);
        
        
        mobility(i)=(L1+L2)/2;
        
        
    end
    for i = 1:size(bedge,1) % loop de fases externas para calcular as mobilidade e fluxo fracional
        v1=bedge(i,1);
        v2=bedge(i,2);
        if bedge(i,5)~=auxflag
            %     lef=bedge(i,3);
            n1=zeros(size(elem,1),1);
            A=zeros(size(elem,1),1);
            n2=zeros(size(elem,1),1);
            A1=zeros(size(elem,1),1);
            for j=esurn2(v1)+1:esurn2(v1+1)
                tt=esurn1(j);
                
                Krw1 = ((S_old(tt) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
                
                Kro1 = ((1 - S_old(tt) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
                
                L11 = Krw1/visc(1) +  Kro1/visc(2); % mobilidade total no elemento esquerdo
                
                A(tt)=elemarea(tt);
                
                n1(tt)=A(tt)*L11;
                
            end
            L1=sum(n1)/sum(A);
            
            for jj=esurn2(v2)+1:esurn2(v2+1)
                tt=esurn1(jj);
                
                Krw1 = ((S_old(tt) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
                
                Kro1 = ((1 - S_old(tt) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
                
                L22 = Krw1/visc(1) +  Kro1/visc(2); % mobilidade total no elemento esquerdo
                
                A1(tt)=elemarea(tt);
                
                n2(tt)=A1(tt)*L22;
                
                
            end
            L2=sum(n2)/sum(A1);
            
            mobility(i+size(inedge,1))=(L1+L2)/2;
        else
            
            Krw1 = ((S_cont - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
            
            Kro1 = ((1 - S_cont - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
            mobility(i+size(inedge,1)) = Krw1/visc(1) +  Kro1/visc(2); % mobilidade total no elemento esquerdo
        end
        
    end
end
end
