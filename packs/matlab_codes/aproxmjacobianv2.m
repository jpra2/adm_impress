function [J]=aproxmjacobianv2(Fk,p_old,nflagface,nflagno,w,s,...
    parameter,weightDMP,kmap,fonte,auxface,mobility,Hesq, Kde, Kn, Kt, ...
    Ded,calnormface,wells,v,M_old,pinterp_new,gravresult,...
        gravrate,gravno,gravelem,gravface,grav_elem_escalar,wg)
global elem
nelem=size(elem,1);
M=full(M_old)~=0;
M=M+0;
J=zeros(nelem,nelem);
x=p_old;
I=eye(nelem);

nlin=size(v,1);

for ielem=1:nlin
    m=0;
    for icol = 1:size(v,2) 
        if v(ielem,icol) ~= 0
            x= x+ 1e-3*sqrt(norm(M_old(:,v(ielem,icol))))*I(:,v(ielem,icol));
            %x= x+ 1e-3*sqrt(norm(norm(p_old)))*I(:,v(ielem,icol)); % Modifiquei aqui!!!!!
            m=m+1;
        end
    end
    
    % Calculo da matriz global
    [auxM,auxRHS]=globalmatrix(x,pinterp_new,0,nflagface,nflagno...
        ,parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
        mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,...
        gravrate,gravno,gravelem,gravface,grav_elem_escalar,wg);
    
    % f(xk+1)
    Fkk= auxM*x - auxRHS;
    
    % Montagem da matriz Jacobiano por método Diferencia Finita
    d(:,ielem) = (Fkk(:) - Fk(:));
    
    for var = 1:m
        if v(ielem,var)~=0
           J(:,v(ielem,var)) = d(:,ielem).*M(:,v(ielem,var))/(1e-3*sqrt(norm(M_old(:,v(ielem,var))))); % modifiquei aqui!!!!
           % J(:,v(ielem,var)) = d(:,ielem).*M(:,v(ielem,var))/(1e-3*sqrt(norm(p_old))); % modifiquei aqui!!!!
        end
    end
    
    % Atualiza o vetor "x"
    x=p_old;
    
end
end