
function [VPI,oilrecovery,cumulateoil,watercut]=reportproduction(S_old,...
    wells,f_elem,cont,VPI,oilrecovery,cumulateoil,watercut,q,d_t)
global elemarea
if max(wells)~=0
    if(cont~=1)
        sumaux1 = 0;
        sumaux2 = 0;
        sumaux3 = 0;
        sumvol  = 0;
        WC=0;
        for iwell = 1:size(wells,1)
            if wells(iwell,2) == 2 % no poço produtor
                sumvol = sumvol + elemarea(wells(iwell,1));
                
                sumaux3 = sumaux3 + f_elem(wells(iwell,1))*elemarea(wells(iwell,1));
                
                f_o = 1 - f_elem(wells(iwell,1));
                
                sumaux2 = sumaux2 + f_o*d_t*elemarea(wells(iwell,1));
                
                WC=WC+S_old(wells(iwell,1))*elemarea(wells(iwell,1));
                
            elseif wells(iwell,2) == 1 % no poço injetor
                
                sumaux1 = sumaux1 + f_elem(wells(iwell,1))*q(wells(iwell,1))*d_t/sum(elemarea);
                
            end
        end
        VPI(cont) = VPI(cont-1) + (sumaux1);
        
        cumulateoil(cont) = cumulateoil(cont-1) + sumaux2/sumvol;
        
        oilrecovery(cont) = 1 - (sumaux3/sumvol);
        
        watercut(cont)= WC/sumvol;
        
    end
else

 oilrecovery=0;
 cumulateoil=0;
 watercut=0;
 VPI=VPI+d_t;
end
end