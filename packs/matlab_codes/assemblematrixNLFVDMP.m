function [M,I]=assemblematrixNLFVDMP(p,pinterp,gamma,nflag,parameter,kmap,...
    fonte,weightDMP,auxface,wells,mobility)

global inedge coord bedge bcflag elem elemarea benchmark

valuemin=1e-16;
I=zeros(size(elem,1),1);
M=zeros(size(elem,1),size(elem,1));
auxmobility1=mobility(1:size(inedge,1),1);
auxmobility2=mobility((size(inedge,1)+1):(size(inedge,1)+size(bedge,1)),1);
mobility(1:size(bedge,1),1)=auxmobility2;
mobility((size(bedge,1)+1):(size(inedge,1)+size(bedge,1)),1)=auxmobility1;
%% fonte
I=I+fonte;

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
            auxmobility=mobility(ifacont);
            %-------------------------------------------------------------%
            % somando 1
            ifacelef1=parameter(1,3,ifacont);
            auxparameter1=parameter(1,1,ifacont);
            [M,I]=tratmentcontourlfvHP(ifacelef1,parameter,nflag,...
                normcont,mobility,auxparameter1,auxmobility,lef,weightDMP,M,I);
            %-------------------------------------------------------------%
            % somando 2
            ifacelef2=parameter(1,4,ifacont);
            auxparameter2=parameter(1,2,ifacont);
            [M,I]=tratmentcontourlfvHP(ifacelef2,parameter,nflag,...
                normcont,mobility,auxparameter2,auxmobility,lef,weightDMP,M,I);
            
            M(lef,lef)=M(lef,lef)+ auxmobility*normcont*(auxparameter1 + auxparameter2);
        end
        
    end
end
x=1000;
y=1000;
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
    
    % Fluxo parcial do elemento a esquerda
    F1=norma*(parameter(1,1,ifactual)*(p(lef)-pinterp(ifacelef1))+ parameter(1,2,ifactual)*(p(lef)-pinterp(ifacelef2)));
    
    % Fluxo parcial do elemento a direita
    F2=norma*(parameter(2,1,ifactual)*(p(rel)-pinterp(ifacerel1))+ parameter(2,2,ifactual)*(p(rel)-pinterp(ifacerel2)));
    
    % caso 1 F1*F2<=0  pag. 2594 Sheng and yuan, 2011
    if F2*F1<0 || F2*F1==0
        % caso 1.1 e caso 1.2, idea adoptado de Gao and Wu, 2013
        mu1=(abs(F2)+valuemin)/(abs(F1)+abs(F2)+2*valuemin);
        
        mu2=(abs(F1)+valuemin)/(abs(F1)+abs(F2)+2*valuemin);
        
        M(lef,lef)=M(lef,lef) + 2*mu1*norma*(parameter(1,1,ifactual)+parameter(1,2,ifactual))*mobility(ifactual);
        M(rel,rel)=M(rel,rel) + 2*mu2*norma*(parameter(2,1,ifactual)+parameter(2,2,ifactual))*mobility(ifactual);
        %somando 1 LEF
        ifacelef1=parameter(1,3,ifactual);
        termo0=2*mu1*norma*parameter(1,1,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMP(M,I,ifacelef1,nflag,lef,parameter,weightDMP,termo0);
        % somando 2 LEF
        
        ifacelef2=parameter(1,4,ifactual);
        termo0=2*mu1*norma*parameter(1,2,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMP(M,I,ifacelef2,nflag,lef,parameter,weightDMP,termo0);
        % somando 1 REL
        
        ifacerel1=parameter(2,3,ifactual);
        termo0=2*mu2*norma*parameter(2,1,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMP(M,I,ifacerel1,nflag,rel,parameter,weightDMP,termo0);
        % somando 2 REL
        
        ifacerel2=parameter(2,4,ifactual);
        termo0=2*mu2*norma*parameter(2,2,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMP(M,I,ifacerel2,nflag,rel,parameter,weightDMP,termo0);
    else
        % Calculo das contribuições do elemento a esquerda
        weightlef=weightDMP(ifactual-size(bedge,1),1);
        weightrel=weightDMP(ifactual-size(bedge,1),2);
        
        % faço esta escolha ja que um das faces pode coincidir com a face
        % em questão
        if ifacelef1==ifactual
            auxparameter1=parameter(1,1,ifactual);
            auxparameter11=parameter(1,2,ifactual);
            auxfaceopost1= ifacelef2;
			% Fluxo parcial do elemento a esquerda
        F1b=norma*(auxparameter11*(p(lef)-pinterp(auxfaceopost1)));
        elseif ifacelef2==ifactual
            
            auxparameter1=parameter(1,2,ifactual);
            auxparameter11=parameter(1,1,ifactual);
            auxfaceopost1= ifacelef1;
			% Fluxo parcial do elemento a esquerda
        F1b=norma*(auxparameter11*(p(lef)-pinterp(auxfaceopost1)));
        else
            % pode ser o caso que a face em questão seja nenhum dos outros
            % faces de interpolação
            auxparameter1=0;
            x=0;
			% Fluxo parcial do elemento a esquerda
        F1b=norma*parameter(1,1,ifactual)*(p(lef)-pinterp(ifacelef1))+ norma*parameter(1,2,ifactual)*(p(lef)-pinterp(ifacelef2));
        end
        %-----------------------------------------------------------------%
        
        % faço esta escolha ja que um das faces pode coincidir com a face
        % em questão
        if ifacerel1==ifactual
            auxparameter2=parameter(2,1,ifactual);
            auxparameter22=parameter(2,2,ifactual);
            auxfaceopost2= ifacerel2;
			% Fluxo parcial do elemento a direita
        F2b=norma*(auxparameter22*(p(rel)-pinterp(auxfaceopost2)));
        elseif ifacerel2==ifactual
            
            auxparameter2=parameter(2,2,ifactual);
            auxparameter22=parameter(2,1,ifactual);
            auxfaceopost2= ifacerel1;
			% Fluxo parcial do elemento a direita
        F2b=norma*(auxparameter22*(p(rel)-pinterp(auxfaceopost2)));
        else
            % pode ser o caso que a face em questão seja nenhum dos outros
            % faces de interpolação
            auxparameter2=0;
            % Fluxo parcial do elemento a direita
        F2b=norma*parameter(2,1,ifactual)*(p(rel)-pinterp(ifacerel1)) + norma*parameter(2,2,ifactual)*(p(rel)-pinterp(ifacerel2));
            y=0;
        end
                
        % calculo dos mu's
        
        mulefb=(abs(F2b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        
        murelb=(abs(F1b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        
        
        
        if (F1b*F2b<0 || F1b*F2b==0) && (x~=0 && y~=0) % caso IV
            alfa=(mulefb*auxparameter1*weightrel + murelb*auxparameter2*weightlef);
            
            % contribuição no lado esquerdo
            M(lef,lef)=M(lef,lef) + mobility(ifactual)*norma*(alfa+2*mulefb*auxparameter11);
            
            M(lef,rel)=M(lef,rel) - mobility(ifactual)*norma*alfa;
            
             % somando 1 LEF
            ifacelef1= auxfaceopost1;
            termo0=2*mulefb*norma*auxparameter11*mobility(ifactual);
          
            [M,I]=auxmatrixassembleDMP(M,I,ifacelef1,nflag,lef,parameter,weightDMP,termo0);
            
            % contribuição no lado direito
            
            M(rel,rel)=M(rel,rel) + mobility(ifactual)*norma*(alfa+2*murelb*auxparameter22);
            
            M(rel,lef)=M(rel,lef) - mobility(ifactual)*norma*alfa;
                
            
            % somando 1 REL
            ifacerel1= auxfaceopost2;
            termo0=2*murelb*norma*auxparameter22*mobility(ifactual);
            
            [M,I]=auxmatrixassembleDMP(M,I,ifacerel1,nflag,rel,parameter,weightDMP,termo0);
			
        elseif (F1b*F2b<0 || F1b*F2b==0) && (x==0 && y~=0) % caso III
		
		         M(lef,lef)=M(lef,lef) + norma*(mulefb*(parameter(1,1,ifactual)+parameter(1,2,ifactual)) + murelb*auxparameter2*weightlef)*mobility(ifactual);
				 M(rel,lef)=M(rel,lef) - norma*(mulefb*(parameter(1,1,ifactual)+parameter(1,2,ifactual)) + murelb*auxparameter2*weightlef)*mobility(ifactual);
				 
                 M(rel,rel)=M(rel,rel) + murelb*norma*(auxparameter22+auxparameter2*weightlef)*mobility(ifactual);
				 M(lef,rel)=M(lef,rel) - murelb*norma*(auxparameter22+auxparameter2*weightlef)*mobility(ifactual);
				 

                 %somando 1 LEF
                 ifacelef1=parameter(1,3,ifactual);
                 termo0=mulefb*norma*parameter(1,1,ifactual)*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef1,nflag,lef,rel,parameter,weightDMP,termo0);
				 
                 % somando 2 LEF
                 ifacelef2=parameter(1,4,ifactual);
                 termo0=mulefb*norma*parameter(1,2,ifactual)*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef2,nflag,lef,rel,parameter,weightDMP,termo0);
			   			              
                 % somando 1 REL
                 ifacerel1= auxfaceopost2;
                 termo0=murelb*norma*auxparameter22*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel1,nflag,rel,lef,parameter,weightDMP,termo0);
			   		
		elseif (F1b*F2b<0 || F1b*F2b==0) && (x~=0 && y==0) % caso II
		
                 % contribuição do elemento a esquerda
		         M(lef,lef)=M(lef,lef) + mulefb*norma*(auxparameter11+auxparameter1*weightrel)*mobility(ifactual);
                 
                 M(lef,rel)=M(lef,rel) - norma*(murelb*(parameter(2,1,ifactual)+parameter(2,2,ifactual)) + mulefb*auxparameter1*weightrel)*mobility(ifactual);
                 
                 % contribuição do elemento a direita
                 M(rel,rel)=M(rel,rel) + norma*(murelb*(parameter(2,1,ifactual)+parameter(2,2,ifactual)) + mulefb*auxparameter1*weightrel)*mobility(ifactual);

                 
                 M(rel,lef)=M(rel,lef) - mulefb*norma*(auxparameter11+auxparameter1*weightrel)*mobility(ifactual);

                 %somando 1 LEF
                 ifacelef1=auxfaceopost1;
                 termo0=mulefb*norma*auxparameter11*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef1,nflag,lef,rel,parameter,weightDMP,termo0);
                 			   			              
                 % somando 1 REL
                 ifacerel1= parameter(2,3,ifactual);
                 termo0=murelb*norma*parameter(2,1,ifactual)*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel1,nflag,rel,lef,parameter,weightDMP,termo0);
		         % somando 2 REL
				 ifacerel2= parameter(2,4,ifactual);
                 termo0=murelb*norma*parameter(2,2,ifactual)*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel2,nflag,rel,lef,parameter,weightDMP,termo0);
		
		
		elseif (F1b*F2b<0 || F1b*F2b==0) && (x==0 && y==0) % caso I
                 % contribuição no elemento a esquerda
                 M(lef,lef)=M(lef,lef) + mulefb*norma*(parameter(1,1,ifactual)+parameter(1,2,ifactual))*mobility(ifactual);
                 M(lef,rel)=M(lef,rel) - murelb*norma*(parameter(2,1,ifactual)+parameter(2,2,ifactual))*mobility(ifactual);
                 
                 % contribuição no elemento a direita
                 M(rel,rel)=M(rel,rel) + murelb*norma*(parameter(2,1,ifactual)+parameter(2,2,ifactual))*mobility(ifactual);
			     M(rel,lef)=M(rel,lef) - mulefb*norma*(parameter(1,1,ifactual)+parameter(1,2,ifactual))*mobility(ifactual);
				 			     
                 %somando 1 LEF
                 ifacelef1=parameter(1,3,ifactual);
                 termo0=mulefb*norma*parameter(1,1,ifactual)*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef1,nflag,lef,rel,parameter,weightDMP,termo0);
                 %somando 2 LEF
                 ifacelef2=parameter(1,4,ifactual);
                 termo0=mulefb*norma*parameter(1,2,ifactual)*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef2,nflag,lef,rel,parameter,weightDMP,termo0);                 
				 
                 % somando 1 REL
                 ifacerel1= parameter(2,3,ifactual);
                 termo0=murelb*norma*parameter(2,1,ifactual)*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel1,nflag,rel,lef,parameter,weightDMP,termo0);
		         % somando 2 REL
				 ifacerel2= parameter(2,4,ifactual);
                 termo0=murelb*norma*parameter(2,2,ifactual)*mobility(ifactual);
                 [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel2,nflag,rel,lef,parameter,weightDMP,termo0);
		
        else
            
            alfa=norma*(mulefb*auxparameter1*weightrel + murelb*auxparameter2*weightlef);
            
            M(lef,lef)=M(lef,lef) + mobility(ifactual)*alfa;
            M(lef,rel)=M(lef,rel) - mobility(ifactual)*alfa;
            
            M(rel,rel)=M(rel,rel) + mobility(ifactual)*alfa;
            M(rel,lef)=M(rel,lef) - mobility(ifactual)*alfa;
            
        end
        
    end
    x=1000;
    y=1000;
    
end
switch benchmark
    case 'gaowu6'
        % malha 23x23
        % M(357,:)=0*M(357,:);
        % M(357,357)=1;
        % I(357)=1;
        % M(173,:)=0*M(173,:);
        % M(173,173)=1;
        % I(173)=0;
        % malha 11x11
        M(83,:)=0*M(83,:);
        M(83,83)=1;
        I(83)=1;
        M(39,:)=0*M(39,:);
        M(39,39)=1;
        I(39)=0;
end
end