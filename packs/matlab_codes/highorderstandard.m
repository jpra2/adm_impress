function [S_old] = highorderstandard(S_old,influx,dt,esuel,wells,totalflux,...
    f_elem,S_cont,Sat_max_min,w,bound,upsilon,kappa,nw,no,auxflag,smetodo)
global elemarea satlimit inedge visc elem bedge pormap  centelem bcflag
% reconstrução dos gradientes utilizando os pesos de minimos quadrados
[g]=grad(S_old,w,esuel,elem,bcflag,bedge,bound);

RHS = zeros(size(S_old,1),1);

for i = 1:size(inedge,1)
    
    v1 = inedge(i,1);
    v2 = inedge(i,2);
    
    outro=inedge(i,4);    %indice do elemento a direita da face i
    aqui=inedge(i,3);    %indice do elemento a esquerda da face i
        
    S_aqui = S_old(aqui);
    S_outro = S_old(outro);
    
    %% calculo dos gradientes centrales
    r_outro = (centelem(aqui,:)-centelem(outro,:))';
    gradS_cent_outro = S_aqui - S_outro; % OK
    
    prod_outro=dot(g(outro,:)',r_outro);
    
    gradS_upwind_outro = 2*prod_outro - gradS_cent_outro; % OK
    
    
    r_aqui = (centelem(outro,:)-centelem(aqui,:))';
    gradS_cent_aqui = S_outro - S_aqui; % OK
    
    prod_aqui=dot(g(aqui,:)',r_aqui);
    
    gradS_upwind_aqui = 2*prod_aqui - gradS_cent_aqui; % OK
   
    %% calculo dos limitadores de Van Albada
    
    limit1= max(0,(2*gradS_upwind_aqui*gradS_cent_aqui+ 1e-16)/(gradS_upwind_aqui^2 + gradS_cent_aqui^2 + 1e-16));
    
    limit2= max(0,(2*gradS_upwind_outro*gradS_cent_outro+ 1e-16)/(gradS_upwind_outro^2 + gradS_cent_outro^2 + 1e-16 ));
    
     %% utilizando o limitador de woodfield ou somente Van Albada
     if strcmp(smetodo,'HOMFV')
         
         [phi_aqui]=limiterwoodf(Sat_max_min, aqui,upsilon, S_old);
         [phi_outro]=limiterwoodf(Sat_max_min, outro,upsilon, S_old);
         
         SD_L = S_aqui  + phi_aqui*limit1*0.25*((1-kappa)*gradS_upwind_aqui+(1+kappa)*gradS_cent_aqui);
         SD_R = S_outro + phi_outro*limit2*0.25*((1-kappa)*gradS_upwind_outro+(1+kappa)*gradS_cent_outro);
     else
         SD_L = S_aqui  + limit1*0.25*((1-kappa)*gradS_upwind_aqui+(1+kappa)*gradS_cent_aqui);
         SD_R = S_outro +  limit2*0.25*((1-kappa)*gradS_upwind_outro+(1+kappa)*gradS_cent_outro);
     end
    %%
    % calculo do fluxo fracional na interface correspondente ao elemento
    % esquerdo
    Krw = ((SD_L - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
    
    Kro = (((1 - SD_L) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
    
    L_L = Krw/visc(1) +  Kro/visc(2);  
    
    ffractional_L = ( Krw/visc(1)) /L_L; 
    % calculo do fluxo fracional na interface conrrespondente ao elemento
    % direita
    Krw = ((SD_R - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
    
    Kro = (((1 - SD_R) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
    
    L_R = Krw/visc(1) +  Kro/visc(2);  
    
    ffractional_R = ( Krw/visc(1)) /L_R; % fluxo fracional da agua no elemento esquerdo
    %% solver de Riemann
    ve_mais  =(influx(i+size(bedge,1)) + abs(influx(i+size(bedge,1))))/2;
    ve_menos =(influx(i+size(bedge,1)) - abs(influx(i+size(bedge,1))))/2;
    RHS(outro) = RHS(outro) + ve_mais*ffractional_L + ve_menos*ffractional_R;
    RHS(aqui)  = RHS(aqui)  - ve_mais*ffractional_L - ve_menos*ffractional_R;
end
%% calculo dos RHS nos poços produtores

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
        if bedge(ifacont,5)==auxflag
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
% calculo das saturações no elemento
for i = 1:size(S_old,1)
    if S_old(i)~=S_cont
        S_old(i,1) = S_old(i,1) + (dt*RHS(i))/(pormap(1)*elemarea(i));
    end
end
end