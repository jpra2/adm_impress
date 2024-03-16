function [M,I]=assemblematrixlfvHPv3(parameter,fonte,nflagface,...
    weightDMP,wells,mobility,gravresult,gravrate,gravelem,gravface)
global inedge coord bedge bcflag elem elemarea gravitational strategy
% incialização das matrizes
I=sparse(size(elem,1),1);
M=sparse(size(elem,1),size(elem,1));
% gravresult:div*flux_g
I=I+fonte;

m=0;
auxmobility=1;
for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef)- normcont*bcflag(r,2);
    else
        if strcmp(gravitational,'yes')
            if strcmp(strategy,'starnoni')
                m=gravrate(ifacont);
            elseif strcmp(strategy,'inhouse')
               m=0;
            end
            
        end
        %-------------------------------------------------------------%
        % somando 1
        ifacelef1=parameter(1,3,ifacont);
        auxparameter1=parameter(1,1,ifacont);
        [M,I]=tratmentcontourlfvHP(ifacelef1,parameter,nflagface,...
            normcont,mobility,auxparameter1,auxmobility,lef,weightDMP,M,I);
        %-------------------------------------------------------------%
        % somando 2
        ifacelef2=parameter(1,4,ifacont);
        auxparameter2=parameter(1,2,ifacont);
        [M,I]=tratmentcontourlfvHP(ifacelef2,parameter,nflagface,...
            normcont,mobility,auxparameter2,auxmobility,lef,weightDMP,M,I);
        
        M(lef,lef)=M(lef,lef)+ auxmobility*normcont*(auxparameter1 + auxparameter2);
        I(lef)=I(lef)+m;
       
    end
    
end

%% Montagem da matriz global

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    % calculo da norma do vetor "vd1"
    norma=norm(vd1);
    % Calculo dos fluxos parciais
    ifactual=iface+size(bedge,1);
    
    %-----------------------------------------------------------------%
    % faces de que contem os pontos de interpolação correspondente a
    % elemento a esquerda
    
    ifacelef1=parameter(1,3,ifactual);
    ifacelef2=parameter(1,4,ifactual);
    % faces de que contem os pontos de interpolação correspondente a
    % elemento a direita
    ifacerel1= parameter(2,3,ifactual);
    ifacerel2= parameter(2,4,ifactual);
    
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
    
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*murel*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    ARR=norma*mulef*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    %     M(lef,lef)=M(lef,lef)+ mobility(ifactual)*ALL;
    %     M(lef,rel)=M(lef,rel)- mobility(ifactual)*ARR;
    %     % contribuição da transmisibilidade no elemento direita
    %     M(rel,rel)=M(rel,rel)+ mobility(ifactual)*ARR;
    %     M(rel,lef)=M(rel,lef)- mobility(ifactual)*ALL;
    
    M(lef,lef)=M(lef,lef)+ 1*ALL;
    M(lef,rel)=M(lef,rel)- 1*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ 1*ARR;
    M(rel,lef)=M(rel,lef)- 1*ALL;
    
    % contribuições do elemento a esquerda
    %------------------------ somando 1 ----------------------------------%
    %termo0=mobility(ifactual)*norma*murel*parameter(1,1,ifactual);
    termo0=1*norma*murel*parameter(1,1,ifactual);
    if ifacelef1<size(bedge,1) || ifacelef1==size(bedge,1)
        if nflagface(ifacelef1,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            pressureface1=nflagface(ifacelef1,2);
            I(lef)=I(lef)+  termo0*pressureface1;
            I(rel)=I(rel)-  termo0*pressureface1;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacelef1"
            % a ideia principal nesta rutina é isolar a pressão na face "ifacelef1"
            % e escrever eem funcção das pressões de face internas ou
            % Dirichlet ou simplesmente em funcação da pressão da VC.
            
            normcont=norm(coord(bedge(ifacelef1,1),:)-coord(bedge(ifacelef1,2),:));
            x=bcflag(:,1)==bedge(ifacelef1,5);
            r=find(x==1);
            
            % calcula o fluxo na face "ifacelef1"
            fluxoN=normcont*bcflag(r,2);
            
            % retorna as faces que formam os eixos auxiliares
            auxifacelef1=parameter(1,3,ifacelef1);
            auxifacelef2=parameter(1,4,ifacelef1);
            
            % retorna as faces os coeficientes correpondente ao face "ifacelef1"
            ksi1=parameter(1,1,ifacelef1);
            ksi2=parameter(1,2,ifacelef1);
            % identifica as faces "atual" que neste caso é "ifacelef1" e a
            % outra face oposto
            if auxifacelef1==ifacelef1
                faceoposto=auxifacelef2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacelef1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            % calcula a mobilidade na face "ifacelef1"
            %mobil=mobility(ifacelef1);
            mobil=1;
            % verifica se a face oposto a "ifacelef1" pertence ao contorno
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    % caso I. Quando a face oposto pertence ao contorno da
                    % malha
                    % atribui a pressão da face de Dirichlet
                    pressure= nflagface(faceoposto,2);
                    % o coeficiente "atualksi" pentence ao face "ifacelef1"
                    termo1=((atualksi+opostoksi)/atualksi);
                    
                    termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    I(lef)=I(lef)-termo0*termo2;
                    % contribuição no elemeto a esquerda
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    I(rel)=I(rel)+termo0*termo2;
                    
                    
                else
                    
                    %mobilO=mobility(faceoposto);
                    mobilO=1;
                    
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    % caso II. Quando a face oposto pertence ao contorno de Neumann
                    
                    if auxifacelef1==faceoposto
                        atualksO=parameter(1,1,faceoposto);
                        opostksO=parameter(1,2,faceoposto);
                    else
                        atualksO=parameter(1,2,faceoposto);
                        opostksO=parameter(1,1,faceoposto);
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    
                    
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(mobilO*normcontopost)) - opostksO*(fluxoN/(mobil*normcont)))/sumatotaldenom;
                    end
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    I(lef)=I(lef)- termo0*termo2;
                    
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    I(rel)=I(rel)+ termo0*termo2;
                    
                end
                
            else
                % caso III. Quando a face oposto pertence ao interior da
                % malha
                termo1=(ksi1+ksi2)/(atualksi);
                
                termo2=(fluxoN/(mobil*normcont*atualksi));
                
                % Calculo das contribuições do elemento a esquerda
                auxweightlef1=weightDMP(faceoposto-size(bedge,1),1);
                auxweightlef2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                if auxlef==lef
                    auxelematual=auxlef;
                    auxelemopost=auxrel;
                    pesatual=auxweightlef1;
                    pesopost=auxweightlef2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweightlef2;
                    pesopost=auxweightlef1;
                    auxelemopost=auxlef;
                end
                
                termo3= (opostoksi/atualksi);
                
                
                % contribuição no elemeto a direita
                M(lef,auxelematual)=M(lef,auxelematual)+ termo0*(termo3*pesatual-termo1);
                
                M(lef,auxelemopost)=M(lef,auxelemopost)+ termo0*termo3*pesopost;
                I(lef)=I(lef)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(rel,auxelematual)=M(rel,auxelematual)- termo0*(termo3*pesatual-termo1);
                
                M(rel,auxelemopost)=M(rel,auxelemopost)- termo0*termo3*pesopost;
                I(rel)=I(rel)+termo0*termo2;
                
            end
            
        end
    else
        % Calculo das contribuições do elemento a esquerda
        auxweightlef1=weightDMP(ifacelef1-size(bedge,1),1);
        auxweightlef2=weightDMP(ifacelef1-size(bedge,1),2);
        
        auxlef=weightDMP(ifacelef1-size(bedge,1),3);
        auxrel=weightDMP(ifacelef1-size(bedge,1),4);
        
        % contribuição do elemento a esquerda
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])- termo0*[auxweightlef1,auxweightlef2];
        
        % contribuição do elemento a direita
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])+ termo0*[auxweightlef1,auxweightlef2];
        
    end
    
    %--------------------------- somando 2 -------------------------------%
    
    %termo0=mobility(ifactual)*norma*murel*parameter(1,2,ifactual);
    termo0=norma*murel*parameter(1,2,ifactual);
    if ifacelef2<size(bedge,1) || ifacelef2==size(bedge,1)
        if nflagface(ifacelef2,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            pressureface2=nflagface(ifacelef2,2);
            I(lef)=I(lef)+ termo0*pressureface2;
            I(rel)=I(rel)- termo0*pressureface2;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacelef2"
            % a ideia principal nesta rutina é isolar a pressão na face "ifacelef2"
            % e escrever eem funcção das pressões de face internas ou
            % Dirichlet ou simplesmente em funcação da pressão da VC.
            
            normcont=norm(coord(bedge(ifacelef2,1),:)-coord(bedge(ifacelef2,2),:));
            x=bcflag(:,1)==bedge(ifacelef2,5);
            r=find(x==1);
            
            % calcula o fluxo na face "ifacelef2"
            fluxoN=normcont*bcflag(r,2);
            
            % retorna as faces que formam os eixos auxiliares
            auxifacelef1=parameter(1,3,ifacelef2);
            auxifacelef2=parameter(1,4,ifacelef2);
            
            % retorna as faces os coeficientes correpondente ao face "ifacelef2"
            ksi1=parameter(1,1,ifacelef2);
            ksi2=parameter(1,2,ifacelef2);
            % identifica as faces "atual" que neste caso é "ifacelef2" e a
            % outra face oposto
            if auxifacelef1==ifacelef2
                faceoposto=auxifacelef2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacelef1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            % calcula a mobilidade na face "ifacelef1"
            mobil=1;
            %mobil=mobility(ifacelef2);
            % verifica se a face oposto a "ifacelef1" pertence ao contorno
            
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    % caso I. Quando a face oposto pertence ao contorno da
                    % malha
                    % atribui a pressão da face de Dirichlet
                    pressure= nflagface(faceoposto,2);
                    % o coeficiente "atualksi" pentence ao face "ifacelef1"
                    termo1=((atualksi+opostoksi)/atualksi);
                    
                    termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    I(lef)=I(lef)-termo0*termo2;
                    % conntribuição no elemento a direita
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    I(rel)=I(rel)+termo0*termo2;
                    
                else
                    %mobilO=mobility(faceoposto);
                    mobilO=1;
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(ifacelef1,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    % caso II. Quando a face oposto pertence ao contorno de Neumann
                    
                    if auxifacelef1==faceoposto
                        atualksO=parameter(1,1,faceoposto);
                        opostksO=parameter(1,2,faceoposto);
                    else
                        atualksO=parameter(1,2,faceoposto);
                        opostksO=parameter(1,1,faceoposto);
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    
                    
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(mobilO*normcontopost)) - opostksO*(fluxoN/(mobil*normcont)))/sumatotaldenom;
                    end
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    I(lef)=I(lef)+ termo0*termo2;
                    
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    I(lef)=I(lef)- termo0*termo2;
                    
                end
                
            else
                % caso III. Quando a face oposto pertence ao interior da
                % malha
                termo1=(ksi1+ksi2)/(atualksi);
                
                termo2=(fluxoN/(mobil*normcont*atualksi));
                
                % Calculo das contribuições do elemento a esquerda
                auxweightlef1=weightDMP(faceoposto-size(bedge,1),1);
                auxweightlef2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                if auxlef==lef
                    auxelematual=auxlef;
                    auxelemopost=auxrel;
                    pesatual=auxweightlef1;
                    pesopost=auxweightlef2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweightlef2;
                    pesopost=auxweightlef1;
                    auxelemopost=auxlef;
                end
                
                termo3= (opostoksi/atualksi);
                
                
                % contribuição no elemeto a direita
                M(lef,auxelematual)=M(lef,auxelematual)+ termo0*(termo3*pesatual-termo1);
                
                M(lef,auxelemopost)=M(lef,auxelemopost)+ termo0*termo3*pesopost;
                I(lef)=I(lef)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(rel,auxelematual)=M(rel,auxelematual)- termo0*(termo3*pesatual-termo1);
                
                M(rel,auxelemopost)=M(rel,auxelemopost)- termo0*termo3*pesopost;
                I(rel)=I(rel)+termo0*termo2;
                
            end
            
        end
    else
        % Calculo das contribuições do elemento a esquerda
        auxweightlef1=weightDMP(ifacelef2-size(bedge,1),1);
        auxweightlef2=weightDMP(ifacelef2-size(bedge,1),2);
        
        auxlef=weightDMP(ifacelef2-size(bedge,1),3);
        auxrel=weightDMP(ifacelef2-size(bedge,1),4);
        
        % contribuição do elemento a esquerda
        M(lef,[auxlef, auxrel])=  M(lef,[auxlef,auxrel])- termo0*[auxweightlef1,auxweightlef2];
        
        % contribuição do elemento a direita
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])+ termo0*[auxweightlef1,auxweightlef2];
        
    end
    %%
    % contribuições do elemento a direita
    % somando 1
    %termo0=mobility(ifactual)*norma*mulef*parameter(2,1,ifactual);
    termo0=norma*mulef*parameter(2,1,ifactual);
    if ifacerel1<size(bedge,1) || ifacerel1==size(bedge,1)
        if nflagface(ifacerel1,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            pressurefacerel1=nflagface(ifacerel1,2);
            I(rel)=I(rel) + termo0*pressurefacerel1;
            I(lef)=I(lef) - termo0*pressurefacerel1;
        else % quando a face é contorno de Neumann
            % agora a face em questão é "ifacerel1"
            % a ideia principal nesta rutina é isolar a pressão na face "ifacerel1"
            % e escrever eem funcção das pressões de face internas ou
            % Dirichlet ou simplesmente em funcação da pressão da VC.
            
            
            normcont=norm(coord(bedge(ifacerel1,1),:)-coord(bedge(ifacerel1,2),:));
            x=bcflag(:,1)==bedge(ifacerel1,5);
            r=find(x==1);
            fluxoN=normcont*bcflag(r,2);
            % lembre-se que a face "ifacerel1" pertence ao contorno, então
            % as faces que definem os eixos auxiliares são
            auxifacerel1=parameter(1,3,ifacerel1);
            auxifacerel2=parameter(1,4,ifacerel1);
            % um deles é mesmo "ifacerel1" e o outro e oposto.
            % as coeficientes que representam a cada face
            ksi1=parameter(1,1,ifacerel1);
            ksi2=parameter(1,2,ifacerel1);
            % aqui calculamos quem é face atual "ifacerel1" e o face oposto
            % e também os coeficientes adequados
            if auxifacerel1==ifacerel1
                faceoposto=auxifacerel2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacerel1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            % calculo da mobilidade da face atual "ifacerel1"
            %mobil=mobility(ifacerel1);
            mobil=1;
            % vejamos a qual face pertence a face oposto
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    % caso I. Quando a face oposto pertence ao contorno da
                    % malha atribui a pressão da face de Dirichlet
                    pressure= nflagface(faceoposto,2);
                    
                    termo1=(atualksi+opostoksi)/atualksi;
                    
                    termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(rel,rel)=M(rel,rel)- termo0*termo1;
                    I(rel)=I(rel)-termo0*termo2;
                    
                    M(lef,rel)=M(lef,rel)+ termo0*termo1;
                    I(lef)=I(lef)+termo0*termo2;
                    
                    
                else
                    %mobilO=mobility(faceoposto);
                    mobilO=1;
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    % caso II. Quando a face oposto pertence ao contorno de Neumann
                    ksi1a=parameter(1,1,faceoposto);
                    ksi2b=parameter(1,2,faceoposto);
                    
                    if auxifacelef1==faceoposto
                        atualksO=ksi1a;
                        opostksO=ksi2b;
                    else
                        atualksO=ksi2b;
                        opostksO=ksi1a;
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    
                    
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(mobilO*normcontopost))...
                     - opostksO*(fluxoN/(mobil*normcont)))/sumatotaldenom;
                        
                    end
                    M(rel,rel)=M(rel,rel)- termo0*termo1;
                    I(rel)=I(rel)+ termo0*termo2;
                    
                    M(lef,rel)=M(lef,rel)+ termo0*termo1;
                    I(lef)=I(lef)- termo0*termo2;
                    
                end
                
            else
                % caso III. Quando a face oposto pertence à interior da
                % malha
                termo1=(opostoksi+atualksi)/(atualksi);
                
                termo2=(fluxoN/(mobil*normcont*atualksi));
                
                termo3= (opostoksi/atualksi);
                
                % lembre-se que a face "ifacerel1" pertence ao contorno.
                % a face oposto a "ifacerel1" é chamado "faceoposto" e
                % pertence ao face interior, então a pressão nessa deve
                % ser interpolado para isso calculamos os pesos e
                % elementos que facem parte da face "oposto"
                
                auxweight1=weightDMP(faceoposto-size(bedge,1),1);
                auxweight2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                
                % lembrando que um "auxlef" ou "auxrel" é "rel"
                if auxlef==rel
                    auxelematual=auxlef;
                    pesatual=auxweight1;
                    
                    auxelemopost=auxrel;
                    pesopost=auxweight2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweight2;
                    
                    pesopost=auxweight1;
                    auxelemopost=auxlef;
                end
                
                % contribuição no elemeto a direita
                M(rel,auxelematual)=M(rel,auxelematual)+ termo0*(termo3*pesatual-termo1 );
                
                M(rel,auxelemopost)=M(rel,auxelemopost)+ termo0*termo3*pesopost;
                
                I(rel)=I(rel)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(lef,auxelematual)=M(lef,auxelematual)- termo0*(termo3*pesatual-termo1 );
                
                M(lef,auxelemopost)=M(lef,auxelemopost)- termo0*termo3*pesopost;
                I(lef)=I(lef)+termo0*termo2;  
            end
            
        end
    else
        % Calculo das contribuições do elemento a esquerda
        auxweightrel1=weightDMP(ifacerel1-size(bedge,1),1);
        auxweightrel2=weightDMP(ifacerel1-size(bedge,1),2);
        
        auxlef=weightDMP(ifacerel1-size(bedge,1),3);
        auxrel=weightDMP(ifacerel1-size(bedge,1),4);

        % contribuição do elemento a esquerda
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])+ termo0*[auxweightrel1,auxweightrel2];
        % contribuição do elemento a direita
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])- termo0*[auxweightrel1,auxweightrel2];
    end
    % somando 2
    %termo0=mobility(ifactual)*norma*mulef*parameter(2,2,ifactual);
    termo0=norma*mulef*parameter(2,2,ifactual);
    if ifacerel2<size(bedge,1) || ifacerel2==size(bedge,1)
        if nflagface(ifacerel2,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            
            pressurefacerel2=nflagface(ifacerel2,2);
            I(rel)=I(rel) + termo0*pressurefacerel2;
            I(lef)=I(lef) - termo0*pressurefacerel2;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacerel1"
            % a ideia principal nesta rutina é isolar a pressão na face "ifacerel2"
            % e escrever eem funcção das pressões de face internas ou
            % Dirichlet ou simplesmente em funcação da pressão da VC.
            
            
            normcont=norm(coord(bedge(ifacerel2,1),:)-coord(bedge(ifacerel2,2),:));
            x=bcflag(:,1)==bedge(ifacerel2,5);
            r=find(x==1);
            fluxoN=normcont*bcflag(r,2);
            % lembre-se que a face "ifacerel1" pertence ao contorno, então
            % as faces que definem os eixos auxiliares são
            auxifacerel1=parameter(1,3,ifacerel2);
            auxifacerel2=parameter(1,4,ifacerel2);
            % um deles é mesmo "ifacerel2" e o outro e oposto.
            % as coeficientes que representam a cada face
            ksi1=parameter(1,1,ifacerel2);
            ksi2=parameter(1,2,ifacerel2);
            % aqui calculamos quem é face atual "ifacerel2" e o face oposto
            % e também os coeficientes adequados
            if auxifacerel1==ifacerel2
                faceoposto=auxifacerel2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacerel1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            % calculo da mobilidade da face atual "ifacerel2"
            %mobil=mobility(ifacerel2);
            mobil=1;
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    % caso I. Quando a face oposto pertence ao contorno da
                    % malha atribui a pressão da face de Dirichlet
                    pressure= nflagface(faceoposto,2);
                    
                    termo1=(atualksi+opostoksi)/atualksi;
                    
                    termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(rel,rel)=M(rel,rel)- termo0*termo1;
                    I(rel)=I(rel)-termo0*termo2;
                    
                    M(lef,rel)=M(lef,rel)+ termo0*termo1;
                    I(lef)=I(lef)+termo0*termo2;
                    
                else
                    mobilO=mobility(faceoposto);
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    % caso II. Quando a face oposto pertence ao contorno de Neumann
                    ksi1a=parameter(1,1,faceoposto);
                    ksi2b=parameter(1,2,faceoposto);
                    
                    if auxifacelef1==faceoposto
                        atualksO=ksi1a;
                        opostksO=ksi2b;
                    else
                        atualksO=ksi2b;
                        opostksO=ksi1a;
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    
                    
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(mobilO*normcontopost))...
                      - opostksO*(fluxoN/(mobil*normcont)))/sumatotaldenom;
                    end
                    M(rel,rel)=M(rel,rel)- termo0*termo1;
                    I(rel)=I(rel)+ termo0*termo2;
                    
                    M(lef,rel)=M(lef,rel)+ termo0*termo1;
                    I(lef)=I(lef)- termo0*termo2;
                    
                end
                
            else
                % caso III. Quando a face oposto pertence à interior da
                % malha
                termo1=(opostoksi+atualksi)/(atualksi);
                
                termo2=(fluxoN/(mobil*normcont*atualksi));
                
                termo3= (opostoksi/atualksi);
                
                % lembre-se que a face "ifacerel2" pertence ao contorno.
                % a face oposto a "ifacerel2" é chamado "faceoposto" e
                % pertence ao face interior, então a pressão nessa deve
                % ser interpolado para isso calculamos os pesos e
                % elementos que facem parte da face "oposto"
                
                auxweight1=weightDMP(faceoposto-size(bedge,1),1);
                auxweight2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                
                % lembrando que um "auxlef" ou "auxrel" é "rel"
                if auxlef==rel
                    auxelematual=auxlef;
                    pesatual=auxweight1;
                    
                    auxelemopost=auxrel;
                    pesopost=auxweight2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweight2;
                    
                    pesopost=auxweight1;
                    auxelemopost=auxlef;
                end
                
                % contribuição no elemeto a direita
                M(rel,auxelematual)=M(rel,auxelematual)+ termo0*(termo3*pesatual-termo1 );
                
                M(rel,auxelemopost)=M(rel,auxelemopost)+ termo0*termo3*pesopost;
                
                I(rel)=I(rel)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(lef,auxelematual)=M(lef,auxelematual)- termo0*(termo3*pesatual-termo1 );
                
                M(lef,auxelemopost)=M(lef,auxelemopost)- termo0*termo3*pesopost;
                I(lef)=I(lef)+termo0*termo2;
                
            end
            
        end
        
    else
        % Calculo das contribuições do elemento a esquerda
        auxweightrel1=weightDMP(ifacerel2-size(bedge,1),1);
        auxweightrel2=weightDMP(ifacerel2-size(bedge,1),2);
        
        auxlef=weightDMP(ifacerel2-size(bedge,1),3);
        auxrel=weightDMP(ifacerel2-size(bedge,1),4);
        
        % contribuição do elemento a direita
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])+ termo0*[auxweightrel1,auxweightrel2];
        
        % contribuição do elemento a esquerda
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])- termo0*[auxweightrel1,auxweightrel2];
        
    end
    if strcmp(gravitational,'yes')
        
        if strcmp(strategy,'starnoni')
            m=gravrate(size(bedge,1)+iface);
        elseif strcmp(strategy,'inhouse')
            m=0;
        end
        
        I(lef)=I(lef)+m;
       
        I(rel)=I(rel)-m;
       
    end
end
%
% % adequação da matriz nos poços produtores
% if max(max(wells))~=0
%     for iw = 1:size(wells,1)
%         if wells(iw,2)==2 %produtor
%             M(wells(iw,1),:)=wells(2,6)*M(wells(iw,1),:);
%             M(wells(iw,1),wells(iw,1))=1;
%             I(wells(iw,1))=0;
%         end
%     end
% end
%% malha 23x23
%M(357,:)=0*M(357,:);
%M(357,357)=1;
%I(357)=1;
%M(173,:)=0*M(173,:);
%M(173,173)=1;
%I(173)=0;
%% malha 11x11
% M(83,:)=0*M(83,:);
% M(83,83)=1;
% I(83)=1;
% M(39,:)=0*M(39,:);
% M(39,39)=1;
% I(39)=0;
end