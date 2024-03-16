function [M,I]=assemblematrixlfvHPv2(w,s,parameter,fonte,nflagface,nflagno,weightDMP,auxface,wells,mobility)
global inedge coord bedge bcflag elem elemarea esurn1 esurn2
% incialização das matrizes
I=zeros(size(elem,1),1);
M=zeros(size(elem,1),size(elem,1));
auxmobility1=mobility(1:size(inedge,1),1);
auxmobility2=mobility((size(inedge,1)+1):(size(inedge,1)+size(bedge,1)),1);
mobility(1:size(bedge,1),1)=auxmobility2;
mobility((size(bedge,1)+1):(size(inedge,1)+size(bedge,1)),1)=auxmobility1;
% fonte
I=I+fonte.*elemarea;
if max(wells)~=0
    sumvol=0;
    for iw = 1:size(wells,1)
        
        if wells(iw,2)==1            % injetor
            I(wells(iw,1))= 1*elemarea(wells(iw,1));        % injeta um m3 de agua por dia (d)
            sumvol=sumvol+ elemarea(wells(iw,1));
        end
    end
    I=I./sumvol;
else
    for ifacont=1:size(bedge,1)
        lef=bedge(ifacont,3);
        
        normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
        
        if bedge(ifacont,5)>200
            x=bcflag(:,1)==bedge(ifacont,5);
            r=find(x==1);
            I(lef)=I(lef)- normcont*bcflag(r,2);
        else
            %-------------------------------------------------------------%
            % somando 1
            if parameter(1,3,ifacont)<size(bedge,1) || ifacelef==size(bedge,1)
                if nflagface(parameter(1,3,ifacont),1)<200 
                    % quando pertence ao contorno de Dirichlet
                    pressureface1=nflagface(parameter(1,3,ifacont),2);
                    Transmicont1= parameter(1,2,ifacont)*pressureface1*normcont;
                    
                else % quando a face é contorno de Neumann
                    
                    no1=bedge(parameter(1,3,ifacont),1);
                    no2=bedge(parameter(1,3,ifacont),2);
                    if nflagno(no1,1)<200
                        
                        I(lef)=I(lef)+0.5*mobility(ifacont)*normcont*parameter(1,1,ifacont)*nflagno(no1,2);
                    else
                        for j=1:(esurn2(no1+1)-esurn2(no1))
                            
                            post_cont=esurn2(no1)+j;
                            
                            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) - 0.5*mobility(ifacont)*normcont*parameter(1,1,ifacont)*w(post_cont);
                            
                        end
                        
                    end
                    
                    if nflagno(no1,1)<200
                        
                        I(lef)=I(lef)+0.5*mobility(ifacont)*normcont*parameter(1,1,ifacont)*nflagno(no2,2);
                        
                    else
                        for j=1:(esurn2(no2+1)-esurn2(no2))
                            
                            post_cont=esurn2(no2)+j;
                            
                            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) - 0.5*mobility(ifacont)*normcont*parameter(1,1,ifacont)*w(post_cont);
                            
                        end
                    end
                    
                    
                end
            else % quando pertence ao face interior
                
                pressureface1=0;
                % implementação da matriz global no contorno
                [ Transmicont1]=auxassemblematrixcontourLINEAR(parameter(1,3,ifacont),...
                    parameter(1,1,ifacont),lef,normcont,weightDMP,bedge);
                
                elemavaliar1=auxface(ifacont,1);
                
                M(lef,elemavaliar1)=M(lef,elemavaliar1)-mobility(ifacont)*Transmicont1;
            end
            %-------------------------------------------------------------%
            % somando 2
            if parameter(1,4,ifacont)<size(bedge,1) || ifacelef==size(bedge,1)
                if nflagface(parameter(1,4,ifacont),1)<200
                    pressureface2=nflagface(parameter(1,4,ifacont),2);
                    
                    Transmicont2= parameter(1,2,ifacont)*pressureface2*normcont;
                    
                else
                    no1=bedge(parameter(1,4,ifacont),1);
                    no2=bedge(parameter(1,4,ifacont),2);
                    if nflagno(no1,1)<200
                        I(lef)=I(lef) + mobility(ifacont)*normcont*parameter(1,2,ifacont)*nflagno(no1,2);
                    else
                        for j=1:(esurn2(no1+1)-esurn2(no1))
                            
                            post_cont=esurn2(no1)+j;
                            
                            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) - 0.5*mobility(ifacont)*normcont*parameter(1,2,ifacont)*w(post_cont);
                            
                        end
                        
                        
                    end
                    
                    if nflagno(no1,1)<200
                        I(lef)=I(lef) + mobility(ifacont)*normcont*parameter(1,2,ifacont)*nflagno(no2,2);
                        
                    else
                        for j=1:(esurn2(no2+1)-esurn2(no2))
                            
                            post_cont=esurn2(no2)+j;
                            
                            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) - 0.5*mobility(ifacont)*normcont*parameter(1,2,ifacont)*w(post_cont);
                            
                        end
                    end
                    
                end
            else
                pressureface2=0;
                [Transmicont2]=auxassemblematrixcontourLINEAR(parameter(1,4,ifacont),...
                    parameter(1,2,ifacont),lef,normcont,weightDMP,bedge);
                
                elemavaliar2=auxface(ifacont,2);
                
                M(lef,elemavaliar2)=M(lef,elemavaliar2)-mobility(ifacont)*Transmicont2;
            end
            %__________________________________________________________________
            
            M(lef,lef)=M(lef,lef)+ mobility(ifacont)*(Transmicont2 + Transmicont1);
            
            I(lef)=I(lef) + mobility(ifacont)*(pressureface1*Transmicont1 + pressureface2*Transmicont2);
        end
    end
end
% Montagem da matriz global

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
    M(lef,lef)=M(lef,lef)+ mobility(ifactual)*ALL;
    M(lef,rel)=M(lef,rel)- mobility(ifactual)*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ mobility(ifactual)*ARR;
    M(rel,lef)=M(rel,lef)- mobility(ifactual)*ALL;
    
    
    % contribuições do elemento a esquerda
    %------------------------ somando 1 ----------------------------------%
    if ifacelef1<size(bedge,1) || ifacelef1==size(bedge,1)
        if nflagface(ifacelef1,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            pressureface1=nflagface(ifacelef1,2);
            I(lef)=I(lef)+ mobility(ifactual)*murel*norma*parameter(1,1,ifactual)*pressureface1;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacelef1"
             termo0=mobility(ifactual)*norma*murel*parameter(1,1,ifactual);
             
             normcont=norm(coord(bedge(ifacelef1,1),:)-coord(bedge(ifacelef1,2),:));
             x=bcflag(:,1)==bedge(ifacelef1,5);
             r=find(x==1);
             fluxoN=normcont*bcflag(r,2);
             
            auxifacelef1=parameter(1,3,ifacelef1);
            auxifacelef2=parameter(1,4,ifacelef1);
            
            ksi1=parameter(1,1,ifacelef1);
            ksi2=parameter(1,2,ifacelef1);
            
            if auxifacelef1==ifacelef1 
                faceoposto=auxifacelef2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacelef1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            mobil=mobility(ifacelef1);
            
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    pressure= nflagface(faceoposto,2);
                    % caso I. Quando a face oposto pertence ao face interior da malha
                    termo1=((ksi1+ksi2)/atualksi);
                    
                    termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    
                    I(lef)=I(lef)+termo2;
                
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
                    
                    sumaksiopost= opostoksi*(atualksO+opostksO);
                    sumaksiatual= opostksO*(atualksi + opostoksi);
                    
                    sumatotalnumer= sumaksiopost+ sumaksiatual;
                    sumatotaldenom= atualksi*opostksO + atualksO*opostoksi;
                    
                   termo1=sumatotalnumer/sumatotaldenom;
                   
                   
                   termo2=(-opostksO*(fluxoN/mobil*normcont) +opostoksi*(fluxOpost/mobilO*normcontopost))/sumatotaldenom;
                   
                   M(lef,lef)=M(lef,lef)- termo0*termo1; 
                   M(rel,lef)=M(rel,lef)+ termo0*termo1;
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
                
                % contribuição no elemeto a esquerda
                
                M(rel,auxelematual)=M(rel,auxelematual)+ termo0*(termo1 - termo3*pesatual);
                
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
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])- mobility(ifactual)*norma*murel*[auxweightlef1,auxweightlef2]*parameter(1,1,ifactual);
                
        % contribuição do elemento a direita
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])+ mobility(ifactual)*norma*murel*[auxweightlef1,auxweightlef2]*parameter(1,1,ifactual);
        
        
    end
    
    %--------------------------- somando 2 -------------------------------%
    
    if ifacelef2<size(bedge,1) || ifacelef2==size(bedge,1)
        if nflagface(ifacelef2,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            pressureface2=nflagface(ifacelef2,2);
            I(lef)=I(lef)+ mobility(ifactual)*murel*norma*parameter(1,2,ifactual)*pressureface2;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacelef1"
             termo0=mobility(ifactual)*norma*murel*parameter(1,2,ifactual);
             
             normcont=norm(coord(bedge(ifacelef2,1),:)-coord(bedge(ifacelef2,2),:));
             x=bcflag(:,1)==bedge(ifacelef2,5);
             r=find(x==1);
             fluxoN=normcont*bcflag(r,2);
             
            auxifacelef1=parameter(1,3,ifacelef2);
            auxifacelef2=parameter(1,4,ifacelef2);
            
            ksi1=parameter(1,1,ifacelef2);
            ksi2=parameter(1,2,ifacelef2);
            
            if auxifacelef1==ifacelef2 
                faceoposto=auxifacelef2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacelef1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            mobil=mobility(ifacelef1);
            
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    pressure= nflagface(faceoposto,2);
                    % caso I. Quando a face oposto pertence ao face interior da malha
                    termo1=((ksi1+ksi2)/atualksi);
                    
                    termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    
                    I(lef)=I(lef)+termo2;
                
                else
                    mobilO=mobility(faceoposto);
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                     x=bcflag(:,1)==bedge(ifacelef1,5);
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
                    
                    sumaksiopost= opostoksi*(atualksO+opostksO);
                    sumaksiatual= opostksO*(atualksi + opostoksi);
                    
                    sumatotalnumer= sumaksiopost+ sumaksiatual;
                    sumatotaldenom= atualksi*opostksO + atualksO*opostoksi;
                    
                   termo1=sumatotalnumer/sumatotaldenom;
                   
                   
                   termo2=(-opostksO*(fluxoN/mobil*normcont) +opostoksi*(fluxOpost/mobilO*normcontopost))/sumatotaldenom;
                   
                   M(lef,lef)=M(lef,lef)- termo0*termo1; 
                   M(rel,lef)=M(rel,lef)+ termo0*termo1;
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
                
                % contribuição no elemeto a esquerda
                
                M(rel,auxelematual)=M(rel,auxelematual)+ termo0*(termo1 - termo3*pesatual);
                
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
        M(lef,[auxlef, auxrel])=  M(lef,[auxlef,auxrel])- mobility(ifactual)*norma*murel*[auxweightlef1,auxweightlef2]*parameter(1,2,ifactual);
                
        % contribuição do elemento a direita
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])+ mobility(ifactual)*norma*murel*[auxweightlef1,auxweightlef2]*parameter(1,2,ifactual);
                
    end
    
    % contribuições do elemento a direita
    % somando 1
    if ifacerel1<size(bedge,1) || ifacerel1==size(bedge,1)
        if nflagface(ifacerel1,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
           
            
            pressurefacerel1=nflagface(ifacerel1,2);
            I(rel)=I(rel) + mobility(ifactual)*norma*mulef*parameter(2,1,ifactual)*pressurefacerel1;
        else % quando a face é contorno de Neumann
            % agora a face em questão é "ifacerel1"
             termo0=mobility(ifactual)*norma*murel*parameter(2,1,ifactual);
             
             normcont=norm(coord(bedge(ifacerel1,1),:)-coord(bedge(ifacerel1,2),:));
             x=bcflag(:,1)==bedge(ifacerel1,5);
             r=find(x==1);
             fluxoN=normcont*bcflag(r,2);
             
            auxifacelef1=parameter(1,3,ifacerel1);
            auxifacelef2=parameter(1,4,ifacerel1);
            
            ksi1=parameter(1,1,ifacerel1);
            ksi2=parameter(1,2,ifacerel1);
            
            if auxifacelef1==ifacerel1 
                faceoposto=auxifacelef2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacelef1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            mobil=mobility(ifacelef1);
            
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    pressure= nflagface(faceoposto,2);
                    % caso I. Quando a face oposto pertence ao face interior da malha
                    termo1=((ksi1+ksi2)/atualksi);
                    
                    termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    
                    I(lef)=I(lef)+termo2;
                
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
                    
                    sumaksiopost= opostoksi*(atualksO+opostksO);
                    sumaksiatual= opostksO*(atualksi + opostoksi);
                    
                    sumatotalnumer= sumaksiopost+ sumaksiatual;
                    sumatotaldenom= atualksi*opostksO + atualksO*opostoksi;
                    
                   termo1=sumatotalnumer/sumatotaldenom;
                   
                   
                   termo2=(-opostksO*(fluxoN/mobil*normcont) +opostoksi*(fluxOpost/mobilO*normcontopost))/sumatotaldenom;
                   
                   M(rel,rel)=M(rel,rel)+ termo0*termo1; 
                   M(lef,rel)=M(lef,rel)- termo0*termo1;
                   I(rel)=I(rel)- termo0*termo2;
                   
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
                M(rel,auxelematual)=M(rel,auxelematual)- termo0*(termo1 - termo3*pesatual);
                
                M(rel,auxelemopost)=M(rel,auxelemopost)+ termo0*termo3*pesopost;
                I(rel)=I(rel)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(lef,auxelematual)=M(lef,auxelematual)+ termo0*(termo1 - termo3*pesatual);
                
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
        
        
        % contribuição do elemento a direita
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])+ mobility(ifactual)*norma*mulef*[auxweightrel1,auxweightrel2]*parameter(2,1,ifactual);
        
        
        % contribuição do elemento a esquerda
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])- mobility(ifactual)*norma*mulef*[auxweightrel1,auxweightrel2]*parameter(2,1,ifactual);
                
        
        
    end
    % somando 2
    if ifacerel2<size(bedge,1) || ifacerel2==size(bedge,1)
        if nflagface(ifacerel2,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            
            pressurefacerel2=nflagface(ifacerel2,2);
            I(rel)=I(rel) + mobility(ifactual)*norma*mulef*parameter(2,2,ifactual)*pressurefacerel2;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacelef1"
             termo0=mobility(ifactual)*norma*murel*parameter(2,2,ifactual);
             
             normcont=norm(coord(bedge(ifacerel2,1),:)-coord(bedge(ifacerel2,2),:));
             x=bcflag(:,1)==bedge(ifacerel2,5);
             r=find(x==1);
             fluxoN=normcont*bcflag(r,2);
             
            auxifacelef1=parameter(1,3,ifacerel2);
            auxifacelef2=parameter(1,4,ifacerel2);
            
            ksi1=parameter(1,1,ifacerel2);
            ksi2=parameter(1,2,ifacerel2);
            
            if auxifacelef1==ifacerel2 
                faceoposto=auxifacelef2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacelef1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            mobil=mobility(ifacerel2);
            
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    pressure= nflagface(faceoposto,2);
                    % caso I. Quando a face oposto pertence ao face interior da malha
                    termo1=((ksi1+ksi2)/atualksi);
                    
                    termo2=(fluxoN/(mobil*normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    
                    I(lef)=I(lef)+termo2;
                
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
                    
                    sumaksiopost= opostoksi*(atualksO+opostksO);
                    sumaksiatual= opostksO*(atualksi + opostoksi);
                    
                    sumatotalnumer= sumaksiopost+ sumaksiatual;
                    sumatotaldenom= atualksi*opostksO + atualksO*opostoksi;
                    
                   termo1=sumatotalnumer/sumatotaldenom;
                   
                   
                   termo2=(-opostksO*(fluxoN/mobil*normcont) +opostoksi*(fluxOpost/mobilO*normcontopost))/sumatotaldenom;
                   
                   M(lef,lef)=M(lef,lef)- termo0*termo1; 
                   M(rel,lef)=M(rel,lef)+ termo0*termo1;
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
                if auxlef==rel
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
                M(rel,auxelematual)=M(rel,auxelematual)- termo0*(termo1 - termo3*pesatual);
                
                M(rel,auxelemopost)=M(rel,auxelemopost)+ termo0*termo3*pesopost;
                I(rel)=I(rel)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(lef,auxelematual)=M(lef,auxelematual)+ termo0*(termo1 - termo3*pesatual);
                
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
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])+ mobility(ifactual)*norma*mulef*[auxweightrel1,auxweightrel2]*parameter(2,2,ifactual);
                
        % contribuição do elemento a esquerda
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])- mobility(ifactual)*norma*mulef*[auxweightrel1,auxweightrel2]*parameter(2,2,ifactual);
                
    end
    
end

% adequação da matriz nos poços produtores
if max(wells)~=0
    for iw = 1:size(wells,1)
        if wells(iw,2)==2 %produtor
            M(wells(iw,1),:)=0*M(wells(iw,1),:);
            M(wells(iw,1),wells(iw,1))=1;
            I(wells(iw,1))=0;
        end
    end
end
%% malha 23x23
% M(357,:)=0*M(357,:);
% M(357,357)=1;
% I(357)=1;
% M(173,:)=0*M(173,:);
% M(173,173)=1;
% I(173)=0;
%% malha 11x11
% M(83,:)=0*M(83,:);
% M(83,83)=1;
% I(83)=1;
% M(39,:)=0*M(39,:);
% M(39,39)=1;
% I(39)=0;
end