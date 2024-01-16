function [gradS_cent_menos_i,gradS_upwind_menos_i,gradS_cent_mas_j,gradS_upwind_mas_j ]= calculgradUpCent(i,j,g,S_old)
global centelem
    rij = (centelem(j,:)-centelem(i,:))';
    
    gradS_cent_menos_i = S_old(j,1) - S_old(i,1); % OK
    
    gradS_upwind_menos_i = 2*dot(g(i,:)',rij) - gradS_cent_menos_i; % OK
    
     gradS_cent_mas_j = S_old(j,1) - S_old(i,1); % OK
    
    gradS_upwind_mas_j = 2*dot(g(j,:)',rij) -  gradS_cent_mas_j; % OK

end