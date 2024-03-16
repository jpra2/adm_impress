function [M,I]=tratmentcontourlfvHP(ifacelef,parameter,nflagface,normcont,...
                         mobility,auxparameter,auxmobility,lef,weightDMP,M,I)
 global bedge bcflag coord
termo0=auxmobility*normcont*auxparameter;
% ifacelef1 por pertencer a face interior, á face Neumann
% ou pode ser a propria face em questão (a face de Dirichlet)
if ifacelef<size(bedge,1) || ifacelef==size(bedge,1)
    if nflagface(ifacelef,1)<200 % é a propia face
        % contorno de Dirichlet
        pressureface1=nflagface(ifacelef,2);
        
        I(lef)=I(lef)+ termo0*pressureface1;
    else % quando a face é contorno de Neumann
        % agora a face em questão é "ifacelef1"
        % a ideia principal nesta rutina é isolar a pressão na face "ifacelef1"
        % e escrever eem funcção das pressões de face internas ou
        % Dirichlet ou simplesmente em funcação da pressão da VC.
        
        normcont=norm(coord(bedge(ifacelef,1),:)-coord(bedge(ifacelef,2),:));
        x=bcflag(:,1)==bedge(ifacelef,5);
        r=find(x==1);
        
        % calcula o fluxo na face "ifacelef1"
        fluxoN=normcont*bcflag(r,2);
        
        % retorna as faces que formam os eixos auxiliares
        auxifacelef1=parameter(1,3,ifacelef);
        auxifacelef2=parameter(1,4,ifacelef);
        
        % retorna as faces os coeficientes correpondente ao face "ifacelef1"
        ksi1=parameter(1,1,ifacelef);
        ksi2=parameter(1,2,ifacelef);
        % identifica as faces "atual" que neste caso é "ifacelef1" e a
        % outra face oposto
        if auxifacelef1==ifacelef
            faceoposto=auxifacelef2;
            atualksi=ksi1;
            opostoksi=ksi2;
        else
            faceoposto=auxifacelef1;
            atualksi=ksi2;
            opostoksi=ksi1;
        end
        % calcula a mobilidade na face "ifacelef1"
        %mobil=mobility(ifacelef);
        mobil=1;
        % verifica se a face oposto a "ifacelef1" pertence ao contorno
        if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
            if nflagface(faceoposto,1)<200
                % caso I. Quando a face oposto pertence ao contorno
                % de Dirichlet
                pressure= nflagface(faceoposto,2);
                % o coeficiente "atualksi" pentence ao face "ifacelef1"
                termo1=((atualksi+opostoksi)/atualksi);
                
                termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                
                % contribuição no elemeto a esquerda
                M(lef,lef)=M(lef,lef)- termo0*termo1;
                I(lef)=I(lef)-termo0*termo2;
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
                I(lef)=I(lef)+ termo0*termo2;
                
                
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
            M(lef,auxelematual)=M(lef,auxelematual)- termo0*(termo1 - termo3*pesatual);
            
            M(lef,auxelemopost)=M(lef,auxelemopost)+ termo0*termo3*pesopost;
            I(lef)=I(lef)-termo0*termo2;
            
        end
        
    end
else
    % caso III. Quando a face oposto pertence ao interior da
    % malha   
    % Calculo das contribuições do elemento a esquerda
    auxweightlef1=weightDMP(ifacelef-size(bedge,1),1);
    auxweightlef2=weightDMP(ifacelef-size(bedge,1),2);
    
    auxlef=weightDMP(ifacelef-size(bedge,1),3);
    auxrel=weightDMP(ifacelef-size(bedge,1),4);
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
        
    % contribuição no elemeto a direita
    M(lef,auxelematual)=M(lef,auxelematual)- termo0*pesatual;
    
    M(lef,auxelemopost)=M(lef,auxelemopost)- termo0*pesopost;    
end

end