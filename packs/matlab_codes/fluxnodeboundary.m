function [Flux1,vetor_face_ordenado]= fluxnodeboundary(inedge,no,influx,w,nsurn2,...
    faces_esurn,k,bflux,ielem)

% do v�rtice em quest�o 'no'
vetor_face_ordenado= faces_esurn(no,1:(nsurn2(no+1)-nsurn2(no)));%vetor de faces
% na vizinhan�a do v�rtice v1

%%
% loop sobre os elementos do vetor_face_ordenado: este vetor ordena
% todas as faces que concorren n� em quest�o em sentido
% antihorario.
j=1;
if length(vetor_face_ordenado)==2
    for iface= vetor_face_ordenado % loop de elementos na vizinhan�a do v�rtice em quest�o
        % t: devolve o n�mero ou flag da face
        if size(inedge,1)<iface % iface pertece ao contorno
            ifacecont= iface - size(inedge,1);
            flux_shapeu(j)= bflux(ifacecont);
        else % iface pertece a malha interior
            if inedge(iface,4)==ielem
                if influx(iface)>0
                    flux=-influx(iface);
                else
                    flux=influx(iface);
                end
            else
                if influx(iface)>0
                    flux=influx(iface);
                else
                    flux=-influx(iface);
                end
            end
            
            if (influx(iface)<0 && flux <0) || (influx(iface)>0 && flux>0)
                flux_shapeu(j)=flux;
            else
                flux_shapeu(j)=0;
            end
            
        end
        j=j+1;
        
    end
    Flux1 =  2*w*sum(flux_shapeu);
else
    for iface= vetor_face_ordenado % loop de elementos na vizinhan�a do v�rtice em quest�o
        % t: devolve o n�mero ou flag da face
        if size(inedge,1)<iface % iface pertece ao contorno
            ifacecont= iface - size(inedge,1);
            flux_shapeu(j)= bflux(ifacecont);
        else % iface pertece a malha interior
            if inedge(iface,4)==ielem
                if influx(iface)>0
                    flux=-influx(iface);
                else
                    flux=influx(iface);
                end
            else
                if influx(iface)>0
                    flux=influx(iface);
                else
                    flux=-influx(iface);
                end
            end
            if (influx(iface)<0 && flux <0) || (influx(iface)>0 && flux>0)
                flux_shapeu(j)=flux;
            else
                flux_shapeu(j)=0;
            end
            
        end
        j=j+1;
        
    end
    Flux1 =  w*sum(flux_shapeu);
end

end