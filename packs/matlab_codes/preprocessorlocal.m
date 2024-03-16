function [pointarmonic,parameter,gamma,p_old,tol,nit,nflagface,...
    nflagno,weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface,calnormface,...
    gravresult,gravrate,weight,contrcontor,wg]=preprocessorlocal(kmap,...
    N,gravelem,gravface)
global elem gravitational pmetodo interpol correction typecorrection
% inicializando as variaveis
nflagno=0;
nflagface=0;
pointarmonic=0;
parameter=0;
auxface=0;
weightDMP=0;
Hesq=0;
Kde=0;
Kn=0;
Kt=0;
Ded=0;
calnormface=0;
gravresult=0;
gravrate=0;
weight=0;
contrcontor=0;
%% dados inicialização métodos dos volumes finitos não linear
gamma=0.0;                     % este parametro esta no intervalo [0,1] pode ser utilizado para o método nao linear MPFA
p_old=1*ones(size(elem,1),1);  % inicializando a presao
tol=1e-10;                      % tolerancia para metodos não lineares
nit=2000;                      % numero de iteracoes de Picard
%% calculo do termo gravitacional
if strcmp(gravitational,'yes')
    
    [gravresult,gravrate]=gravitation(kmap,gravelem,gravface);
    
end
%% Calculo dos pesos

if strcmp(pmetodo,'nlfvLPEW')|| strcmp(pmetodo,'nlfvLPS') || ...
        strcmp(pmetodo,'nlfvHP')||strcmp(pmetodo,'nlfvDMPSY')||...
        strcmp(pmetodo,'lfvHP')|| strcmp(pmetodo,'lfvLPEW')|| ...
        strcmp(pmetodo,'mpfad')
    % adequação dos flags de contorno
    nflagface= contflagface;
    % calculo dos pesos que correspondem aos metodos de interpolacao
    if strcmp(interpol,'LPEW1')
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_1(kmap,N);
    elseif strcmp(interpol,'eLPEW2')
        % interpolaca LPEW2 modificado por proposto por Miao e Wu 2021
        [weight,contrcontor] = Pre_ELPEW_2(kmap,N,gravrate);
    elseif strcmp(interpol,'LS')
        [ weight,contrcontor] = LS(kmap);
    elseif strcmp(interpol,'eLS')
        disp('>> falta implementar!')
    else
        % interpolaca LPEW2 proposto por Gao e Wu 2010
        [weight,contrcontor,wg] = Pre_LPEW_2(kmap,N,gravrate,gravelem);
    end
    
end
%% calculo dos variavei inherentes ao metodo
if strcmp(pmetodo,'nlfvLPEW')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter,calnormface]=coefficientLPSangle(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno;
    
elseif strcmp(pmetodo,'nlfvLPS')|| strcmp(pmetodo,'lfvLPEW')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter,calnormface]=coefficientLPSangle(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno;
    if strcmp(pmetodo,'lfvLPEW')
        [weightDMP]=weightnlfvDMP(kmap);
        % outra maneira de calcular os pesos proposto no artigo
        %[weightDMP]=weightlfv(parameter);
    end
    
elseif strcmp(pmetodo,'interpfree')
    [parameter]=coeffinterpfree(kmap,F);
    
elseif strcmp(pmetodo,'nlfvPPS')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter,calnormface]=coefficientLPSangle(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno;
    
    if strcmp(correction,'yes')
        if strcmp(typecorrection,'firstcorrection')
            % correcao utilizando express. simplificada
            [pointarmonic,weightDMP,raioaux]=firstcorrectharmonic(kmap,N);
        elseif strcmp(typecorrection,'secondcorrection')
            % correcao ao ponto medio
            [pointarmonic,weightDMP,raioaux]=secondcorrectharmonic(kmap,N);
        else
            % correcao dos pontos harmonicos segundo Kobaise
            [pointarmonic,weightDMP,raioaux]=thirdcorrectharmonic(kmap);
        end
        
    else
        % calculoa dos pontos harmonicos sem correcao 'express. original'
        [pointarmonic,weightDMP,raioaux]=nocorrectharmonic(kmap);
    end
elseif strcmp(pmetodo,'nlfvDMPSY')|| strcmp(pmetodo,'lfvHP') || ...
        strcmp(pmetodo,'nlfvDMPV1')|| strcmp(pmetodo,'nlfvHP')
    %% faces alrededor de um elemento
    [facelement]=element_face;
    
    if strcmp(correction,'yes')
        if strcmp(typecorrection,'firstcorrection')
            % correcao utilizando express. simplificada
            [pointarmonic,weightDMP,raioaux]=firstcorrectharmonic(kmap,N);
        elseif strcmp(typecorrection,'secondcorrection')
            % correcao ao ponto medio
            [pointarmonic,weightDMP,raioaux]=secondcorrectharmonic(kmap,N);
        else
            % correcao dos pontos harmonicos segundo Kobaise
            [pointarmonic,weightDMP,raioaux]=thirdcorrectharmonic(kmap);
        end
        
    else
        % calculoa dos pontos harmonicos sem correcao 'express. original'
        [pointarmonic,weightDMP,raioaux]=nocorrectharmonic(kmap);
    end
    %% calculo dos parametros ou constantes (ksi)
    % temos usado este parametro durante muito tempo em muitos testes
    [parameter,auxface]=coefficientPPSharmonicpoint(facelement,pointarmonic,kmap,raioaux);
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPSusingHP(kmap,facelement,pointarmonic); %para lfvHP
    % adequação dos flag de face de contorno
    nflagface= contflagface;
    % adequação dos nos flags de contorno
    nflagno= contflagno;
    %calculo de parametros
    if strcmp(pmetodo,'nlfvHP')
        [Hesq, Kde, Kn, Kt, Ded] = Kde_Ded_Kt_Kn( kmap);
    end
elseif strcmp(pmetodo,'mpfad')
    
    % calculo das constantes fisicos-geometrico
    [Hesq, Kde, Kn, Kt, Ded] = Kde_Ded_Kt_Kn(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno;
else
    % calculo das constantes fisicos-geometrico para o TPFA
    [Hesq, Kde, Kn, Kt, Ded] = Kde_Ded_Kt_Kn( kmap);
    
    % adequação dos flags de contorno
    nflagface= contflagface;
end

end