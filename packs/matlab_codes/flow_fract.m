%funçao que calcula: permeabilidades relativas - Krw e Kro
%                    mobilidades - lambda_w e lambda_o
%                    fluxo fracional - f

function [fM] = flow_fract(SM,nw,no)
global satlimit visc
   % calculo das permeabilidade usando modelo Brooks e Corey
    
    % calculando as permeabilidades relativas e a mobilidade do elemento a esquerdo
    
    Krw = ((SM - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
    
    Kro = (((1 - SM) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
    
    L1 = Krw/visc(1) +  Kro/visc(2);  
    
    fM = ( Krw/visc(1)) /L1; % fluxo fracional da agua no elemento esquerdo
end
