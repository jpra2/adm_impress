function [S_old]= firstorderstandard(S_old,influx,totalflux,f_elem,dt,wells,S_cont,auxflag,nw,no)

global elem inedge bedge elemarea pormap satlimit visc

RHS=zeros(size(elem,1),1);

%% aproximação upwind
for iface = 1 : size(inedge,1)
    
    lef = inedge(iface,3); % elemento a esquerda
    rel = inedge(iface,4); % elemento a direita
    ve_mais = (influx(iface+size(bedge,1)) + abs(influx(iface+size(bedge,1))))/2;
    ve_menos = (influx(iface+size(bedge,1)) - abs(influx(iface+size(bedge,1))))/2;
    
    RHS(rel) = RHS(rel) + ve_mais*f_elem(lef) + ve_menos*f_elem(rel);
    RHS(lef)  = RHS(lef)  - ve_mais*f_elem(lef) - ve_menos*f_elem(rel);
    
end
%% calculo dos RHS nos poços produtores
% calculo dos contribuições nos poços
if max(wells)~=0
    for iw=1:size(wells,1)
        
        if wells(iw,2)==2
           
            RHS(wells(iw,1)) = RHS(wells(iw,1)) + f_elem(wells(iw,1))*totalflux(wells(iw,1));
        end
    end
    % calculo dos contribuições nos contornos onde tem fluxo prescrito
    % diferente de zero
else
    for ifacont = 1:size(bedge,1)
        lef=bedge(ifacont,3);
        if bedge(ifacont,5)==auxflag || bedge(ifacont,5)==102
            % calculo do fluxo fracional no contorno com saturação
            % prescrito
            Krw2 = ((S_cont - satlimit(1))/(1-satlimit(1)-satlimit(2)))^nw;
            
            Kro2 = ((1 - S_cont - satlimit(1))/(1-satlimit(1)-satlimit(2)))^no;
            
            L2 = Krw2/visc(1) + Kro2/visc(2);
            
            f_cont = (Krw2/visc(1))/L2;
            
            RHS(lef) = RHS(lef) - f_cont * influx(ifacont);
        elseif bedge(ifacont,5)<200
            
            RHS(bedge(ifacont,3)) = RHS(bedge(ifacont,3)) - f_elem(bedge(ifacont,3))*influx(ifacont);
        end
    end
end
%% calculo da saturação
for i = 1:size(elem,1)
    por=pormap(1);
    if S_old(i,1)~=S_cont
        S_old(i,1) = S_old(i,1) + (dt*RHS(i))/(por*elemarea(i));
    end
end
end