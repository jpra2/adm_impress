function [pressurinterp]=pressureinterp(p,nflagface,nflagno,w,s,...
    parameter,weightDMP,mobility,gravresult)
global inedge coord esurn1 esurn2 bedge bcflag pmetodo

if strcmp(pmetodo,'nlfvDMPSY')
    pressurinterp=zeros(size(inedge,1)+size(bedge,1),1);
    %% interpolação das pressões no pontos armonicos internos
    pressurinterp(size(bedge,1)+1:size(bedge,1)+size(inedge,1),1)=weightDMP(:,1).*p(weightDMP(:,3))+ weightDMP(:,2).*p(weightDMP(:,4));
    %% interpolação das pressões no pontos medios do contorno
    for ifacont=1:size(bedge,1)
        
        norma=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
        
        % Quando o ponto harmonico oposto pertece a face interior da malha ou
        % contorno de Dirichlet
        
        if bedge(ifacont,5)==201
            no1=bedge(ifacont,1);
            no2=bedge(ifacont,2);
            
            nec1=esurn2(no1+1)-esurn2(no1);
            p1=0;
            % calculo da pressão no vértice v1
            for j=1:nec1
                element1=esurn1(esurn2(no1)+j);
                p1=p1+w(esurn2(no1)+j)*p(element1);
            end
            %% pressão no vértice v1
            p1=p1+ s(no1,1);
            % calculo da pressão no vértice v2
            nec2=esurn2(no2+1)- esurn2(no2);
            p2=0;
            for j=1:nec2
                element2=esurn1(esurn2(no2)+j);
                p2=p2+w(esurn2(no2)+j)*p(element2);
            end
            %% pressão no vértice v2
            p2=p2+ s(no2,1);
            
            %% pressão no ponto medio da face
            pressuremedio=(p1+p2)*0.5;
            
            pressurinterp(ifacont,1)=pressuremedio;
        elseif bedge(ifacont,5)==205
            
            A=logical(parameter(1,3,ifacont)==ifacont);
            B=logical(parameter(1,4,ifacont)==ifacont);
            
            auxface=A*parameter(1,4,ifacont) + B*parameter(1,3,ifacont);
            r=A*1+B*2;
            auxr=A*2+B*1;
            
            if auxface>size(bedge,1)
                % pressão oposto na face interna
                pressurinterp(ifacont,1)=  (((parameter(1,2,ifacont)+ parameter(1,1,ifacont))/parameter(1,r,ifacont))*p(bedge(ifacont,3))-...
                    (parameter(1,auxr,ifacont)/parameter(1,r,ifacont))*pressurinterp(auxface) -  (nflagface(ifacont,2)/(norma*parameter(1,r,ifacont))));
            elseif nflagface(auxface,1)<200
                
                % pressão oposto na face de Dirichlet
                pressurinterp(ifacont,1)= (((parameter(1,2,ifacont)+ parameter(1,1,ifacont))/parameter(1,r,ifacont))*p(bedge(ifacont,3))-...
                    (parameter(1,auxr,ifacont)/parameter(1,r,ifacont))*nflagface(auxface,2)-(nflagface(ifacont,2)/(norma*parameter(1,r,ifacont))));
            end
            
        else
            pressurinterp(ifacont,1)= nflagface(ifacont,2);
        end
        
    end
    
elseif strcmp(pmetodo,'nlfvHP')
    
    
    pressurinterp=zeros(size(inedge,1)+size(bedge,1),1);
    %% interpolação das pressões no pontos armonicos internos
    pressurinterp(size(bedge,1)+1:size(bedge,1)+size(inedge,1),1)=weightDMP(:,1).*p(weightDMP(:,3))+ weightDMP(:,2).*p(weightDMP(:,4));
    %% interpolação das pressões no pontos medios do contorno
    for ifacont=1:size(bedge,1)
        
        norma=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
        
       
        % contorno de Neumann
        if bedge(ifacont,5)==201
            no1=bedge(ifacont,1);
            no2=bedge(ifacont,2);
            %% pressao no vértice v1
            % calculo da pressao no vertice v1
            nec1=esurn2(no1+1)-esurn2(no1);
            p1=0;
            
            for j=1:nec1
                element1=esurn1(esurn2(no1)+j);
                p1=p1+w(esurn2(no1)+j)*p(element1);
            end
            
            p1=p1+ s(no1,1);
            %% pressao no vertice v2
            % calculo da pressao no vertice v2
            nec2=esurn2(no2+1)- esurn2(no2);
            p2=0;
            for j=1:nec2
                element2=esurn1(esurn2(no2)+j);
                p2=p2+w(esurn2(no2)+j)*p(element2);
            end
            
            p2=p2+ s(no2,1);
            
            %% pressao no ponto medio da face
            pressuremedio=(p1+p2)*0.5;
            
            pressurinterp(ifacont,1)=pressuremedio;
        elseif bedge(ifacont,5)==205
             % Quando o ponto harmonico oposto pertece a face interior da malha ou
             % contorno de Dirichlet
            A=logical(parameter(1,3,ifacont)==ifacont);
            B=logical(parameter(1,4,ifacont)==ifacont);
            
            auxface=A*parameter(1,4,ifacont) + B*parameter(1,3,ifacont);
            r=A*1+B*2;
            auxr=A*2+B*1;
            
            if auxface>size(bedge,1)
                % pressao oposto na face interna
                pressurinterp(ifacont,1)=  (((parameter(1,2,ifacont)+ parameter(1,1,ifacont))/parameter(1,r,ifacont))*p(bedge(ifacont,3))-...
                    (parameter(1,auxr,ifacont)/parameter(1,r,ifacont))*pressurinterp(auxface) -  (nflagface(ifacont,2)/(norma*parameter(1,r,ifacont))));
            elseif nflagface(auxface,1)<200
                
                % pressao oposto na face de Dirichlet
                pressurinterp(ifacont,1)= (((parameter(1,2,ifacont)+ parameter(1,1,ifacont))/parameter(1,r,ifacont))*p(bedge(ifacont,3))-...
                    (parameter(1,auxr,ifacont)/parameter(1,r,ifacont))*nflagface(auxface,2)-(nflagface(ifacont,2)/(norma*parameter(1,r,ifacont))));
            end
            
        else
            % contorno de Dirichlet
            pressurinterp(ifacont,1)= nflagface(ifacont,2);
        end
        
    end
    
elseif strcmp(pmetodo,'lfvHP')||strcmp(pmetodo,'nlfvDMPV1')
    mobil=1;
    mobilO=1;
    pressurinterp=zeros(size(inedge,1)+size(bedge,1),1);
    %% interpolação das pressões no pontos armonicos internos
    pressurinterp(size(bedge,1)+1:size(bedge,1)+size(inedge,1),1)=weightDMP(:,1).*p(weightDMP(:,3))+ weightDMP(:,2).*p(weightDMP(:,4));
    %% interpolação das pressões no pontos medios do contorno
    for ifacont=1:size(bedge,1)
        lef=bedge(ifacont,3);
        
        % Quando o ponto harmonico oposto pertece a face interior da malha ou
        % contorno de Dirichlet
        
        if bedge(ifacont,5)>200
            auxfacelef1=parameter(1,3,ifacont);
            auxfacelef2=parameter(1,4,ifacont);
            if auxfacelef1==ifacont
                faceoposto=auxfacelef2;
                atualksi=parameter(1,1,ifacont);
                opostoksi=parameter(1,2,ifacont);
            else
                faceoposto=auxfacelef1;
                atualksi=parameter(1,2,ifacont);
                opostoksi=parameter(1,1,ifacont);
            end
            normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
            x=bcflag(:,1)==bedge(ifacont,5);
            r=find(x==1);
            
            % calcula o fluxo na face "ifacelef2"
            fluxoN=normcont*bcflag(r,2);
            %mobil=mobility(ifacont);
            % "faceoposto" é aquele face oposto ao "ifacont"
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                % faceoposto pode ser contorno de Dirichlet
                if nflagface(faceoposto,1)<200
                    %  caso I. Quando a face oposto pertence ao contorno de Dirichlet
                    pressure=nflagface(faceoposto,2);
                    
                    termo1=(opostoksi+atualksi)/atualksi;
                    termo2= fluxoN/(mobil*normcont*atualksi);
                    termo3= opostoksi/atualksi;
                    
                    pressurinterp(ifacont,1)= termo1*p(lef)-termo2-termo3*pressure;
                else
                    %  caso II. Quando a face oposto pertence ao contorno de Neumann
                    %mobilO=mobility(faceoposto);
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    
                    if auxfacelef1==faceoposto
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
                    pressurinterp(ifacont,1)= termo1*p(lef)+termo2;
                    
                end
            else
                % caso III. Quando a face oposto pertence ao interior da
                % malha
                termo1=(opostoksi+atualksi)/atualksi;
                
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
                
                pfaceopost= pesatual*p(auxelematual)+pesopost*p(auxelemopost);
                
                pressurinterp(ifacont,1)= termo1*p(auxelematual)-termo2-termo3*pfaceopost;
            end
        else
            pressurinterp(ifacont,1)= nflagface(ifacont,2);
        end
        
    end
    
elseif strcmp(pmetodo,'nlfvLPEW') || strcmp(pmetodo,'lfvLPEW') || ...
        strcmp(pmetodo,'nlfvPPS') ||strcmp(pmetodo,'mpfad') ||strcmp(pmetodo,'nlfvLPS')
    %% interpolação das pressões nos nós
    for no=1:size(coord,1)
        
        nec1=esurn2(no+1)-esurn2(no);
        p1=0;
        
        if nflagno(no,1) >200
            if no<=length(s)
                for j=1:nec1
                    element1=esurn1(esurn2(no)+j);
                    p1=p1+w(esurn2(no)+j)*p(element1);
                end
                p1=p1+s(no,1);
            else
                for j=1:nec1
                    element1=esurn1(esurn2(no)+j);
                    p1=p1+w(esurn2(no)+j)*p(element1);
                end
            end
            
        else
            p1=nflagno(no,2);
        end
        
        pressurinterp(no,1)=p1;
        
    end
else
    
    pressurinterp=0;
    
end

end