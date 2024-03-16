function [M,I]=auxmatrixassembleDMP(M,I,iface,nflagface,ielem,...
                                     parameter,weightDMP,termo0)
global bedge coord bcflag

if iface<size(bedge,1) || iface==size(bedge,1)
    if nflagface(iface,1)<200
        % neste caso automaticamente introduzimos o valor dada no
        % contorno de Dirichlet
        % neste caso interpolamos pelo lpew os vertices da face
        pressureface1=nflagface(iface,2);
        I(ielem)=I(ielem)+  termo0*pressureface1;
      
    else% quando a face � contorno de Neumann
        % agora a face em quest�o � "ifacelef1"
        % a ideia principal nesta rutina � isolar a press�o na face "ifacelef1"
        % e escrever eem func��o das press�es de face internas ou
        % Dirichlet ou simplesmente em funca��o da press�o da VC.
        
        normcont=norm(coord(bedge(iface,1),:)-coord(bedge(iface,2),:));
        x=bcflag(:,1)==bedge(iface,5);
        r=find(x==1);
        
        % calcula o fluxo na face "ifacelef1"
        fluxoN=normcont*bcflag(r,2);
        
        % retorna as faces que formam os eixos auxiliares
        auxifacelef1=parameter(1,3,iface);
        auxifacelef2=parameter(1,4,iface);
        
        % retorna as faces os coeficientes correpondente ao face "ifacelef1"
        ksi1=parameter(1,1,iface);
        ksi2=parameter(1,2,iface);
        % identifica as faces "atual" que neste caso � "ifacelef1" e a
        % outra face oposto
        if auxifacelef1==iface
            faceoposto=auxifacelef2;
            atualksi=ksi1;
            opostoksi=ksi2;
        else
            faceoposto=auxifacelef1;
            atualksi=ksi2;
            opostoksi=ksi1;
        end
        % calcula a mobilidade na face "ifacelef1"
        mobil=mobility(iface);
        % verifica se a face oposto a "ifacelef1" pertence ao contorno
        if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
            if nflagface(faceoposto,1)<200
                % caso I. Quando a face oposto pertence ao contorno da
                % malha
                % atribui a press�o da face de Dirichlet
                pressure= nflagface(faceoposto,2);
                % o coeficiente "atualksi" pentence ao face "ifacelef1"
                termo1=((atualksi+opostoksi)/atualksi);
                
                termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                
                % contribui��o no elemeto a esquerda
                M(ielem,ielem)=M(ielem,ielem)- termo0*termo1;
                I(ielem)=I(ielem)-termo0*termo2;
          
            else
                
                mobilO=mobility(faceoposto);
                
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
                M(ielem,ielem)=M(ielem,ielem)- termo0*termo1;
                I(ielem)=I(ielem)- termo0*termo2;
                                
            end
            
        else
            % caso III. Quando a face oposto pertence ao interior da
            % malha
            termo1=(ksi1+ksi2)/(atualksi);
            
            termo2=(fluxoN/(mobil*normcont*atualksi));
            
            % Calculo das contribui��es do elemento a esquerda
            auxweightlef1=weightDMP(faceoposto-size(bedge,1),1);
            auxweightlef2=weightDMP(faceoposto-size(bedge,1),2);
            
            auxlef=weightDMP(faceoposto-size(bedge,1),3);
            auxrel=weightDMP(faceoposto-size(bedge,1),4);
            if auxlef==ielem
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
            
            
            % contribui��o no elemeto a direita
            M(ielem,auxelematual)=M(ielem,auxelematual)+ termo0*(termo3*pesatual-termo1);
            
            M(ielem,auxelemopost)=M(ielem,auxelemopost)+ termo0*termo3*pesopost;
            I(ielem)=I(ielem)-termo0*termo2;
            
            
        end
        
    end
else
    % Calculo das contribui��es do elemento a esquerda
    auxweightlef1=weightDMP(iface-size(bedge,1),1);
    auxweightlef2=weightDMP(iface-size(bedge,1),2);
    
    auxlef=weightDMP(iface-size(bedge,1),3);
    auxrel=weightDMP(iface-size(bedge,1),4);
    
    % contribui��o do elemento a esquerda
    M(ielem,[auxlef,auxrel])=  M(ielem,[auxlef,auxrel])- termo0*[auxweightlef1,auxweightlef2];
       
end
end