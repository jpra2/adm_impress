function [M,I]=assemblematrixlfvHP(w,s,parameter,fonte,nflagface,nflagno,weightDMP,auxface,wells,mobility)
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
if max(max(wells))~=0
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
                    pressureface1=nflagface(parameter(1,3,ifacont),2);
                    Transmicont1= parameter(1,2,ifacont)*pressureface1*normcont;
                    
                else
                    % quando a face é contorno de Neumann
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
            else
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
        else
            
            % quando a face é contorno de Neumann
            no1=bedge(ifacelef1,1);
            no2=bedge(ifacelef1,2);
            if nflagno(no1,1)<200
                
                pressureno1=nflagno(no1,2);
                I(lef)=I(lef)+ 0.5*mobility(ifactual)*norma*murel*parameter(1,1,ifactual)*pressureno1;
                
            else
                auxvetor=1:(esurn2(no1+1)-esurn2(no1));
                auxvetor=auxvetor+esurn2(no1);
                auxesurn1=esurn1(auxvetor);
                
                M(lef, auxesurn1)=M(lef,auxesurn1) - 0.5*mobility(ifactual)*norma*murel*parameter(1,1,ifactual)*w(auxvetor);
                
                M(rel, auxesurn1)=M(rel,auxesurn1) + 0.5*mobility(ifactual)*norma*murel*parameter(1,1,ifactual)*w(auxvetor);
                
            end
            
            if nflagno(no2,1)<200
                pressureno2=nflagno(no2,2);
                I(lef)=I(lef)+ 0.5*mobility(ifactual)*norma*murel*parameter(1,1,ifactual)*pressureno2;
                
            else
                
                auxvetor=1:(esurn2(no2+1)-esurn2(no2));
                auxvetor=auxvetor+esurn2(no2);
                auxesurn1=esurn1(auxvetor);

                M(lef, auxesurn1)=M(lef,auxesurn1) - 0.5*mobility(ifactual)*norma*murel*parameter(1,1,ifactual)*w(auxvetor);
                
                M(rel, auxesurn1)=M(rel,auxesurn1) + 0.5*mobility(ifactual)*norma*murel*parameter(1,1,ifactual)*w(auxvetor);

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
        else
            % quando a face é contorno de Neumann
            no1=bedge(ifacelef2,1);
            no2=bedge(ifacelef2,2);
            
            if nflagno(no1,1)<200
                pressureno1=nflagno(no1,2);
                I(lef)=I(lef)+ 0.5*mobility(ifactual)*norma*murel*parameter(1,2,ifactual)*pressureno1;
                
            else
                auxvetor=1:(esurn2(no1+1)-esurn2(no1));
                auxvetor=auxvetor+esurn2(no1);
                auxesurn1=esurn1(auxvetor);
                
                M(lef, auxesurn1)=M(lef,auxesurn1) - 0.5*mobility(ifactual)*norma*murel*parameter(1,2,ifactual)*w(auxvetor);
                
                M(rel, auxesurn1)=M(rel,auxesurn1) + 0.5*mobility(ifactual)*norma*murel*parameter(1,2,ifactual)*w(auxvetor);
                
            end
            if nflagno(no2,1)<200
                pressureno2=nflagno(no2,2);
                
                I(lef)=I(lef)+ 0.5*mobility(ifactual)*norma*murel*parameter(1,2,ifactual)*pressureno2;
                
            else
                auxvetor=1:(esurn2(no2+1)-esurn2(no2));
                auxvetor=auxvetor+esurn2(no2);
                auxesurn1=esurn1(auxvetor);
                
                M(lef, auxesurn1)=M(lef,auxesurn1) - 0.5*mobility(ifactual)*norma*murel*parameter(1,2,ifactual)*w(auxvetor);
                
                M(rel, auxesurn1)=M(rel,auxesurn1) + 0.5*mobility(ifactual)*norma*murel*parameter(1,2,ifactual)*w(auxvetor);
                
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
            % neste caso interpolamos pelo lpew os vertices da face
            
            pressurefacerel1=nflagface(ifacerel1,2);
            I(rel)=I(rel) + mobility(ifactual)*norma*mulef*parameter(2,1,ifactual)*pressurefacerel1;
        else
            
            % quando a face é contorno de Neumann
            no1=bedge(ifacerel1,1);
            no2=bedge(ifacerel1,2);
            
            if nflagno(no1,1)<200
                pressureno1=nflagno(no1,2);
                
                I(rel)=I(rel)+0.5*mobility(ifactual)*norma*mulef*parameter(2,1,ifactual)* pressureno1;
                
            else
                auxvetor=1:(esurn2(no1+1)-esurn2(no1));
                auxvetor=auxvetor+esurn2(no1);
                auxesurn1=esurn1(auxvetor);
                
                M(rel, auxesurn1)=M(rel, auxesurn1) - 0.5*mobility(ifactual)*norma*mulef*parameter(2,1,ifactual)*w(auxvetor);
                
                M(lef, auxesurn1)=M(lef, auxesurn1) + 0.5*mobility(ifactual)*norma*mulef*parameter(2,1,ifactual)*w(auxvetor);
   
            end
            if nflagno(no2,1)<200
                
                pressureno2=nflagno(no2,2);
                
                I(rel)=I(rel)+0.5*mobility(ifactual)*norma*mulef*parameter(2,1,ifactual)* pressureno2;

            else
                auxvetor=1:(esurn2(no2+1)-esurn2(no2));
                auxvetor=auxvetor+esurn2(no2);
                auxesurn1=esurn1(auxvetor);
               
                M(rel, auxesurn1)=M(rel,auxesurn1) - 0.5*mobility(ifactual)*norma*mulef*parameter(2,1,ifactual)*w(auxvetor);
                
                M(lef, auxesurn1)=M(lef,auxesurn1) + 0.5*mobility(ifactual)*norma*mulef*parameter(2,1,ifactual)*w(auxvetor);
                
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
        else
            
            % quando a face é contorno de Neumann
            no1=bedge(ifacerel2,1);
            no2=bedge(ifacerel2,2);
            if nflagno(no1,1)<200
                pressureno1=nflagno(no1,2);
                I(rel)=I(rel)+ 0.5*mobility(ifactual)*mulef*parameter(2,2,ifactual)*norm*pressureno1;
            else
                auxvetor=1:(esurn2(no1+1)-esurn2(no1));
                auxvetor=auxvetor+esurn2(no1);
                auxesurn1=esurn1(auxvetor);
                
                M(rel, auxesurn1)=M(rel,auxesurn1) - 0.5*mobility(ifactual)*mulef*parameter(2,2,ifactual)*norma*w(auxvetor);
                
                M(lef, auxesurn1)=M(rel,auxesurn1) + 0.5*mobility(ifactual)*mulef*parameter(2,2,ifactual)*norma*w(auxvetor);
                
            end
            if nflagno(no2,1)<200
                pressureno2=nflagno(no2,2);
                I(rel)=I(rel)+ 0.5*mobility(ifactual)*mulef*parameter(2,2,ifactual)*norma*pressureno2;
            else
                auxvetor=1:(esurn2(no2+1)-esurn2(no2));
                auxvetor=auxvetor+esurn2(no2);
                auxesurn1=esurn1(auxvetor);
                
                M(rel,auxesurn1)=M(rel,auxesurn1) - 0.5*mobility(ifactual)*mulef*parameter(2,2,ifactual)*norma*w(auxvetor);
                
                M(lef,auxesurn1)=M(lef,auxesurn1) + 0.5*mobility(ifactual)*mulef*parameter(2,2,ifactual)*norma*w(auxvetor);
                
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
if max(max(wells))~=0
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