function [pressure,errorelativo,flowrate,flowresult,tabletol,coercividade]=...
    solverpressure(kmap,nflagface,nflagno,fonte,...
    tol, nit,p_old,mobility,gamma,wells,parameter,...
    Hesq, Kde, Kn, Kt, Ded,weightDMP,auxface,...
    calnormface,gravresult,gravrate,w,s,gravno,gravelem,gravface,grav_elem_escalar,wg)
global iteration pmetodo elemarea strategy
errorelativo=0;
tabletol=0;
coercividade=0;

%% calculo da pressao 
switch pmetodo
    % Calculo da pressao utilizando metodos nao-lineares
    case {'nlfvLPEW', 'nlfvDMPSY','nlfvDMPV1','nlfvHP', 'nlfvPPS','nlfvLPS'}
        % intepolacao dos pontos auxiliarea com o antigo campo de pressao
        [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,parameter,...
            weightDMP,mobility);
        % calculo da matrizes globlais iniciais
        [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
            parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
            mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
        if strcmp(iteration,'AA')
            % calculo das variavel pressao 
            tic
            [pressure,tabletol,iter,ciclos]=picardAA(M_old,RHS_old,nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
            toc
        elseif strcmp(iteration,'fullpicard')
            % calculo das variavel pressao 
            tic
            [pressure,tabletol,iter,ciclos]=fullpicard(M_old,RHS_old,...
                nit,tol,kmap,parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
            toc
        elseif strcmp(iteration,'iterbroyden')
            % calculo das variavel pressao 
             tic                 
            [pressure, iter,ciclos,tolerancia]=broyden(M_old,RHS_old,...
                p_old,tol,kmap,parameter,w,s,nflagface,fonte,gamma,nflagno,...
                weightDMP,auxface,calnormface,wells,mobility,gravresult,...
                gravrate,gravno,gravelem,gravface,grav_elem_escalar,wg,pinterp);
            toc
        elseif strcmp(iteration,'RRE')
            % calculo das variavel pressao 
            tic
            [pressure,tabletol,iter,ciclos]=picardRRE(M_old,RHS_old,nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
            toc
        elseif strcmp(iteration,'MPE')
            % calculo das variavel pressao 
            tic
            [pressure,tabletol,iter,ciclos]=picardMPE(M_old,RHS_old,nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
            toc
        elseif strcmp(iteration, 'iterdiscretnewton')
            
            p_old1=M_old\RHS_old;
            % interpolação nos nós ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,nflagno,w,s,...
                parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            % resolvedor de pressão pelo método de Newton-Discreto
            [pressure,iter,ciclos,tolerancia]=iterdiscretnewton(M_old1,...
                RHS_old1,M_old,RHS_old,nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,...
                weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt,...
                Ded,calnormface,p_old1);
            
            
        elseif strcmp(iteration, 'iterhybrid')
            
            p_old1=M_old\RHS_old;
            
            % interpolação nos nós ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,w,s,parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,w,s,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            % solver pressure pelo método hybrido
            [pressure,iter,ciclos,tolerancia]=iterhybrid(M_old1,RHS_old1,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,p_old1,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded);
            
        elseif strcmp(iteration, 'JFNK')
            
            [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,...
                parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            p_old1=M_old\RHS_old;
            % calculo do residuo
            R0=norm(M_old*p_old-RHS_old);
            
            % interpolação nos nós ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,nflagno,w,s,...
                parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,...
                nflagface,nflagno,...
                parameter,kmap,fonte,w,s,weightDMP,...
                auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            % calculo da pressão
            [pressure,iter,ciclos,tolerancia]= JFNK1(tol,kmap,parameter,w,...
                s,nflagface,fonte,gamma,...
                nflagno,M_old1,RHS_old1,p_old1,R0,weightDMP,...
                auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
        elseif strcmp(iteration, 'fsolver')
            tic
           [pressure,tabletol,iter,ciclos]=fsolver(M_old,RHS_old,...
                nit,tol,kmap,parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
           toc
        end
        
        %  % informacoes da simulacao
        niteracoes=iter*ciclos;
        name = pmetodo;
        X = sprintf('>> Pressure field calculation according to the method: %s ',name);
        disp(X)
        if strcmp(iteration,'iterbroyden')
            x=['>> Tolerance:',num2str(tolerancia)];
            disp(x);
        else
            x=['>> Tolerance:',num2str(tabletol(size(tabletol,1),2))];
            disp(x);
        end
        name1=iteration;
        z = sprintf('>> was used the iterative method: %s ',name1);
        disp(z)
        y=['>> Iteration numbers:',num2str(niteracoes)];
        disp(y);
        
    % Calculo da pressao utilizando metodos lineares (MPFAs)    
    case {'lfvHP','lfvLPEW','mpfad','tpfa'}
        % incializando variaveis
        
        pinterp=0;
        % calculo da matriz globlal inicial
        [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
            parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
            mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface,grav_elem_escalar,wg);
        
            if strcmp(strategy,'inhouse1')
                pressure=M_old\RHS_old - elemarea.*grav_elem_escalar(:);
            else
            pressure=M_old\RHS_old;
           end
        % informacoes da simulacao
        tabletol=0;
        name = pmetodo;
        X = sprintf('>> Pressure field calculation according to the method: %s ',name);
        disp(X)  
end
%% calculo das vazoes 
% intepolacao dos pontos auxiliare a com o novo campo de pressao
pinterp=pressureinterp(pressure,nflagface,nflagno,w,s,parameter,weightDMP);
% flowrate: vazoes
% flowresult: balanco de massa em cada elemento
if strcmp(pmetodo,'nlfvDMPSY')
    %implementação do fluxo NLFV-DMP
    [flowrate,flowresult]=flowrateNLFVDMP(pressure, pinterp, parameter,...
        nflagface,kmap,gamma,weightDMP,mobility);
elseif strcmp(pmetodo,'nlfvHP')
    [flowrate,flowresult]=flowrateNLFVHP(pressure, pinterp, parameter,gravrate);
    coercividade=0;
elseif strcmp(pmetodo,'nlfvLPEW')
    %implementação do fluxo NLFV
    [flowrate,flowresult,coercividade]=flowrateNLFV(pressure, pinterp,...
        parameter,mobility,gravrate);
elseif strcmp(pmetodo, 'nlfvLPS')
    %implementação do fluxo NLFV
    [flowrate,flowresult,coercividade]=flowrateNLFVPP(pressure, pinterp,...
        parameter,mobility);
    
elseif strcmp(pmetodo, 'nlfvPPS')
    %implementação do fluxo NLFV
    [flowrate,flowresult,coercividade]=flowrateNLFVPP(pressure, pinterp,...
        parameter,mobility);
    
elseif strcmp(pmetodo, 'lfvHP')
    [flowrate,flowresult]=flowratelfvHP(parameter,weightDMP,mobility,...
        pinterp,pressure,gravrate);
    
elseif strcmp(pmetodo, 'lfvLPEW')
    %calculo das vazões
    [flowrate,flowresult]=flowratelfvLPEW(parameter,weightDMP,mobility,...
        pinterp,pressure);
    
elseif  strcmp(pmetodo, 'tpfa')
    [flowrate, flowresult]=flowrateTPFA(pressure,Kde,Kn,Hesq,nflagface,...
        mobility,gravresult,gravrate,pinterp);
else
    %calculo das vazões
    [flowrate,flowresult]=calflowrateMPFAD(pressure,w,s,Kde,Ded,Kn,Kt,...
        Hesq,nflagno,1,gravresult,gravrate,pinterp,gravno,gravelem,grav_elem_escalar);
end


end