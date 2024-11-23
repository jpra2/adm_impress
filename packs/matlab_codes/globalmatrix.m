% funcao que faz a montagem da matrizes globais para todos os metodos
% estudado
function [M,I]=globalmatrix(p,pinterp,gamma,nflagface,nflagno,parameter,kmap,...
    fonte,w,s,weightDMP,auxface,wells,mobility,Hesq,...
    Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,gravno,...
    gravelem,gravface,grav_elem_escalar,wg)

            
global pmetodo

if strcmp(pmetodo,'nlfvDMP1')
    
    [M,I]=assemblematrixDMP(pinterp,gamma,nflagface,parameter,kmap,fonte);
elseif strcmp(pmetodo,'nlfvDMP2')
    
    [M,I]=assemblematrixDMPv1(pinterp,gamma,nflagface,parameter,kmap,fonte);
elseif strcmp(pmetodo,'nlfvLPS') || strcmp(pmetodo,'nlfvPPS')
    
    [M,I]=assemblematrixLPSPPS(p,pinterp,parameter,fonte);
elseif strcmp(pmetodo,'nlfvLPEW')
    
    [M,I]=assemblematrixGYZS(pinterp,parameter,fonte,wells,calnormface,gravrate);
    
elseif strcmp(pmetodo,'nlfvHP')
    [M,I]=assemblematrixNLFVHP(pinterp,parameter,fonte,wells,Hesq,...
        Kn,Kt,nflagno,gravrate,gravresult,nflagface);
elseif strcmp(pmetodo,'nlfvDMPSY')
    
    [M,I]=assemblematrixDMPSY(p,pinterp,gamma,nflagface,parameter,kmap,fonte,...
        weightDMP,auxface,wells,mobility,gravrate);
elseif strcmp(pmetodo,'nlfvDMPV1')
    
    [M,I]=assemblematrixNLFVDMP(p,pinterp,gamma,nflagface,parameter,kmap,...
    fonte,weightDMP,auxface,wells,mobility);
elseif strcmp(pmetodo,'lfvLPEW')
    
    [M,I]=assemblematrixlfvLPEW(parameter,fonte,w,s,nflagno,weightDMP,...
        wells,mobility);
elseif strcmp(pmetodo,'lfvHP')
    
    [M,I]=assemblematrixlfvHPv3(parameter,fonte,nflagface,weightDMP,wells,...
        mobility,gravresult,gravrate,gravelem,gravface);
elseif strcmp(pmetodo,'mpfad')
    [ M, I ] = globalmatrixmpfad(w,s, Kde, Ded, Kn, Kt, nflagno, Hesq,...
        fonte,gravresult,gravrate,gravno,gravelem,grav_elem_escalar,wg);
elseif strcmp(pmetodo,'tpfa')
    [ M, I ] = globalmatrixtpfa( Kde, Kn, nflagface, Hesq,gravresult,...
        gravrate,gravno,gravelem,fonte);
   
end
end