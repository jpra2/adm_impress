%funçao que calcula: permeabilidades relativas - Krw e Kro
%                    mobilidades - lambda_w e lambda_o
%                    fluxo fracional - f

function [f_elem] = fractionalflow(S_old,nw,no)
global visc satlimit elem
f_elem=zeros(size(elem,1),1);

% loop de fases internas para calcular as mobilidade e fluxo fracional
for ielem = 1:size(elem,1) 
    % calculo das permeabilidade usando modelo Brooks e Corey
     
    % calculando as permeabilidades relativas e a mobilidade do elemento a esquerdo

    Krw2 = ((S_old(ielem) - satlimit(1))/(1-satlimit(1)-satlimit(2)))^nw;
    
    Kro2 = ((1 - S_old(ielem) - satlimit(1))/(1-satlimit(1)-satlimit(2)))^no;
    
    L2 = Krw2/visc(1) + Kro2/visc(2);
         
    f_elem(ielem) = (Krw2/visc(1))/L2; 
 
end
       
end
