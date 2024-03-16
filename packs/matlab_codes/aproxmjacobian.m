function [J]=aproxmjacobian(Fk,p_old,nflagface,nflagno,w,s,...
    parameter,weightDMP,kmap,fonte,auxface,mobility,Hesq, Kde, Kn, Kt, ...
    Ded,calnormface,wells,pinterp_new,gravresult,gravrate,gravno,...
    gravelem,gravface,grav_elem_escalar,wg)

global elem
nelem=size(elem,1);
J=sparse(nelem,nelem);
x=p_old;
I=eye(nelem);

delta = 1e-3*sqrt(norm(p_old));

for ielem=1:nelem
    %  xi+h
    x = x + delta*I(:,ielem);
    
    % Calculo da matriz global
    [auxM,auxRHS]=globalmatrix(x,pinterp_new,0,nflagface,nflagno...
        ,parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
        mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,...
        gravrate,gravno,gravelem,gravface,grav_elem_escalar,wg);
    % f(xk+1)
    Fkk= auxM*x - auxRHS;
    
    % Montagem da matriz Jacobiano por método Diferencia Finita
    J(1:nelem,ielem)=(Fkk(:)-Fk(:))./(delta);
     %o errotá aq      
         
    % Atualiza o vetor "x"
    x=p_old;
    delta=1e-3*sqrt(norm(x));
end
end