function [p,step,errorelativo,flowrate,flowresult,tabletol]=fullpicard_teste(M_old,RHS_old,nitpicard,tolpicard,kmap,...
    parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,...
    weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface)

%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);

%% inicializando dados para iteração Picard
step=0;
er=1;
contador=0;
while (tolpicard<er || tolpicard==er) && (step<nitpicard)
    %% atualiza iterações
    step=step+1
    
    for i=1:8
        %p_new=M_old\RHS_old;  % inversão sem pivotamento
        [L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-9));
        %[L,U] = lu(M_old);
        %[p_new]=bicgstab(M_old,RHS_old,1e-8,1000,L,U);
        [p_new]=gmres(M_old,RHS_old,10,1e-9,1000,L,U);
        %% plotagem no visit
        %     S=ones(size(p_new,1),1);
        %     postprocessor(p_new,S,step)
        %     p_max=max(p_new)
        %     p_min=min(p_new)
        %% Interpolação das pressões na arestas (faces)
        [pinterp_new]=pressureinterp(p_new,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);
        
        %% Calculo da matriz global
        [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflagface,nflagno...
            ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
        
        %% calculo do erro
        R = norm(M_new*p_new - RHS_old);
        
        if (R0 ~= 0.0)
            er = abs(R/R0)
            if er<tolpicard
                break
            end
        else
            er = 0.0; %exact
        end
        M_old=M_new;
    end
        %% Interpolação das pressões na arestas (faces)
        [pinterp_new]=pressureinterp(p_new,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);
        
        %% Calculo da matriz global
        [M_old,RHS_old]=globalmatrix(p_new,pinterp_new,gamma,nflagface,nflagno...
            ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
         R0=norm(M_old*p_new - RHS_old);
         
    errorelativo(step)=er;
    %M_old=M_new;
    %RHS_old=RHS_new;
    %if rem(step,5)==0 || step==1
    tabletol(contador+1,1:2)=[contador, er];
    contador=contador+1;
    %end
end
p=p_new;
%p=M_old\RHS_old;
pinterp=pressureinterp(p,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);

if strcmp(metodoP,'nlfvDMPSY')
    % implementação do fluxo NLFV-DMP
    [flowrate,flowresult]=flowrateNLFVDMP(p, pinterp, parameter,nflagface,kmap,gamma,weightDMP,mobility);
else
    % implementação do fluxo NLFV
    [flowrate,flowresult]=flowrateNLFV(p, pinterp, parameter,mobility);
end
residuo=er;
niteracoes=step;

name = metodoP;
X = sprintf('Calculo o campo de pressão pelo método: %s ',name);
disp(X)

x=['Erro:',num2str(residuo)];
disp(x);
y=['Número de iterações:',num2str(niteracoes)];
disp(y);

end