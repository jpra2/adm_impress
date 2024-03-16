function [M,I]=assemblematrixDMPGAOWU(p,pinterp,gamma,nflag,parameter,kmap,...
    fonte,benchmark,weightDMP,auxface,wells,mobility)
global inedge coord bedge bcflag elem elemarea


I=zeros(size(elem,1),1);
M=zeros(size(elem,1),size(elem,1));

auxmobility1=mobility(1:size(inedge,1),1);
auxmobility2=mobility((size(inedge,1)+1):(size(inedge,1)+size(bedge,1)),1);
mobility(1:size(bedge,1),1)=auxmobility2;
mobility((size(bedge,1)+1):(size(inedge,1)+size(bedge,1)),1)=auxmobility1;
%% fonte
I=I+fonte.*elemarea;
%%
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
%% Montagem da matriz global
x=1000;
y=1000;
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    % calculod dos pesos de interpolação
    [weightlef,weightrel]=weightnlfvDMP(kmap,iface);
    %% Calculo dos fluxos parciais
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
    weightlef=weightDMP(ifactual-size(bedge,1),1);
    weightrel=weightDMP(ifactual-size(bedge,1),2);
    
    % faço esta escolha ja que um das faces pode coincidir com a face
    % em questão
    if ifacelef1==ifactual
        auxparameter1=parameter(1,1,ifactual);
        auxparameter11=parameter(1,2,ifactual);
        auxfaceopost1= ifacelef2;
        % calculo do fluxo a aesquerda
        partialfluxlef=norma*((1-gamma)*parameter(1,1,ifactual)*(p(lef)-pinterp(ifactual))+ ...
            parameter(1,2,ifactual)*(p(lef)-pinterp(auxfaceopost1)));
        
    elseif ifacelef2==ifactual
        
        auxparameter1=parameter(1,2,ifactual);
        auxparameter11=parameter(1,1,ifactual);
        auxfaceopost1= ifacelef1;
        % calculo do fluxo a aesquerda
        partialfluxlef=norma*((1-gamma)*parameter(1,2,ifactual)*(p(lef)-pinterp(ifactual))+ ...
            parameter(1,1,ifactual)*(p(lef)-pinterp(auxfaceopost1)));
    else
        % pode ser o caso que a face em questão seja nenhum dos outros
        % faces de interpolação
        auxparameter1=0;
        x=0;
        % Fluxo parcial do elemento a esquerda
        partialfluxlef=norma*(parameter(1,1,ifactual)*(p(lef)-pinterp(parameter(1,3,ifactual)))+ ...
            parameter(1,2,ifactual)*(p(lef)-pinterp(parameter(1,4,ifactual))));
    end
    %-----------------------------------------------------------------%
    
    % faço esta escolha ja que um das faces pode coincidir com a face
    % em questão
    if ifacerel1==ifactual
        auxparameter2=parameter(2,1,ifactual);
        auxparameter22=parameter(2,2,ifactual);
        auxfaceopost2= ifacerel2;
        % calculo do fluxo parcial do elemento a direita
        
        partialfluxrel=norma*((1-gamma)*parameter(2,1,ifactual)*(p(rel)-pinterp(ifactual))+ ...
            parameter(2,2,ifactual)*(p(rel)-pinterp(auxparameter22)));
    elseif ifacerel2==ifactual
        
        auxparameter2=parameter(2,2,ifactual);
        auxparameter22=parameter(2,1,ifactual);
        auxfaceopost2= ifacerel1;
        % calculo do fluxo parcial do elemento a direita
        partialfluxrel=norma*((1-gamma)*parameter(2,2,ifactual)*(p(rel)-pinterp(ifactual))+ ...
            parameter(2,1,ifactual)*(p(rel)-pinterp(auxparameter22)));
    else
        % pode ser o caso que a face em questão seja nenhum dos outros
        % faces de interpolação
        auxparameter2=0;
        % Fluxo parcial do elemento a direita
        y=0;
        % calculo do fluxo parcial do elemento a direita
        partialfluxrel=norma*(parameter(2,1,ifactual)*(p(rel)-pinterp(parameter(2,3,ifactual)))+ ...
            parameter(2,2,ifactual)*(p(rel)-pinterp(parameter(2,4,ifactual))));
        
    end
    %% Calculo das constantes da não linearidade (eq. 13) do artigo Gao e Wu, (2013)
    % veja a REMARK 3.2 pag. 313 do mesmo artigo
    
    mulef=abs(partialfluxrel+e-16)/(abs(partialfluxlef)+abs(partialfluxrel) + 2*e-16);
    
    murel=abs(partialfluxlef+e-16)/(abs(partialfluxlef)+abs(partialfluxrel) + 2*e-16);
    
    %% implementação do  beta
    beta=(1-sign(partialfluxlef*partialfluxrel));
    if  x~=0 && y~=0 % caso IV
        
        % calculo das transmisibilidade
        alfa=(1-gamma)*norma*(mulef*auxparameter1*weightrel + murel*auxparameter2*weightlef);
        
        % Contribuição do elemento a esquerda
        M(lef,lef)=M(lef,lef)+mobility(ifactual)*norma*( alfa+(1-gamma)*mulef*beta*(auxparameter1 + auxparameter11 ));
        
        M(lef,rel)=M(lef,rel)-mobility(ifactual)*norma*alfa;
        
        % Contribuição do elemento a direita
        M(rel,rel)=M(rel,rel)+mobility(ifactual)*norma*( alfa+(1-gamma)*murel*beta*(auxparameter2 + auxparameter22 ));
        
        M(rel,lef)=M(rel,lef)-mobility(ifactual)*norma*alfa;
        
        % estas contribuições surgem apartir do fluxo parcial
        % termo 1 LEF
        % "auxparameter1" face atual corerspondente a elemento a esquerda
        termo0=(1-gamma)*mulef*beta*auxparameter1;
        ifacelef1=auxparameter1;
        [M,I]=auxmatrixassembleDMP(M,I,ifacelef1,nflagface,lef,parameter,weightDMP,termo0);
        
        % termo 2 LEF
        % "auxparameter11" face oposta a face atual correspondente a elemento a esquerda
        termo0=mulef*beta*auxparameter11;
        ifacelef2=auxparameter11;
        [M,I]=auxmatrixassembleDMP(M,I,ifacelef2,nflagface,lef,parameter,weightDMP,termo0);
        
        % termo 1 REL
        % "auxparameter2" face atual correspondente a elemento a direita
        termo0=(1-gamma)*murel*beta*auxparameter2;
        ifacerel1=auxparameter2;
        [M,I]=auxmatrixassembleDMP(M,I,ifacerel1,nflagface,rel,parameter,weightDMP,termo0);
        
        % termo 2 REL
        % "auxparameter2" face oposto a face atual correspondente a elemento direita
        
        termo0=murel*beta*auxparameter22;
        ifacerel2=auxparameter22;
        [M,I]=auxmatrixassembleDMP(M,I,ifacerel2,nflagface,rel,parameter,weightDMP,termo0);
        
        % y~=0 --> o elemento a direita recebe a contribuição do elemento a esquerda diretamente
        % x==0 --> o elemento a esquerda não recebe a contribuição do elemento a direita
        
    elseif x==0 && y~=0 % caso III
        
        M(lef,lef)=M(lef,lef) + norma*(mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual)) + murelb*auxparameter2*weightlef)*mobility(ifactual);
        M(rel,lef)=M(rel,lef) - norma*(mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual)) + murelb*auxparameter2*weightlef)*mobility(ifactual);
        
        M(rel,rel)=M(rel,rel) + murel*norma*(auxparameter22+auxparameter2*weightlef)*mobility(ifactual);
        M(lef,rel)=M(lef,rel) - murel*norma*(auxparameter22+auxparameter2*weightlef)*mobility(ifactual);
        
        
        %somando 1 LEF
        ifacelef1=parameter(1,3,ifactual);
        termo0=mulef*norma*parameter(1,1,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef1,nflag,lef,rel,parameter,weightDMP,termo0);
        
        % somando 2 LEF
        ifacelef2=parameter(1,4,ifactual);
        termo0=mulef*norma*parameter(1,2,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef2,nflag,lef,rel,parameter,weightDMP,termo0);
        
        % somando 1 REL
        ifacerel1= auxfaceopost2;
        termo0=murel*norma*auxparameter22*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel1,nflag,rel,lef,parameter,weightDMP,termo0);
        
    elseif x~=0 && y==0 % caso II
        
        M(lef,lef)=M(lef,lef) + mulef*norma*(auxparameter11+auxparameter1*weightrel)*mobility(ifactual);
        M(rel,lef)=M(rel,lef) - mulef*norma*(auxparameter11+auxparameter1*weightrel)*mobility(ifactual);
        
        M(rel,rel)=M(rel,rel) + norma*(murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual)) + mulef*auxparameter1*weightrel)*mobility(ifactual);
        M(lef,rel)=M(lef,rel) - norma*(murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual)) + mulef*auxparameter1*weightrel)*mobility(ifactual);
        %somando 1 LEF
        ifacelef1=auxfaceopost1;
        termo0=mulef*norma*auxparameter11*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef1,nflag,lef,rel,parameter,weightDMP,termo0);
        
        % somando 1 REL
        ifacerel1= parameter(2,3,ifactual);
        termo0=murel*norma*parameter(2,1,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel1,nflag,rel,lef,parameter,weightDMP,termo0);
        % somando 2 REL
        ifacerel2= parameter(2,4,ifactual);
        termo0=murel*norma*parameter(2,2,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel2,nflag,rel,lef,parameter,weightDMP,termo0);
        
        
    elseif x==0 && y==0 % caso I
	    % contribuição no elemento a esquerda
        M(lef,lef)=M(lef,lef) + mulef*norma*(parameter(1,1,ifactual)+parameter(1,2,ifactual))*mobility(ifactual);
		M(lef,rel)=M(lef,rel) - murel*norma*(parameter(2,1,ifactual)+parameter(2,2,ifactual))*mobility(ifactual);
		
		% contribuição no elemento a direita
		M(rel,rel)=M(rel,rel) + murel*norma*(parameter(2,1,ifactual)+parameter(2,2,ifactual))*mobility(ifactual);
        M(rel,lef)=M(rel,lef) - mulef*norma*(parameter(1,1,ifactual)+parameter(1,2,ifactual))*mobility(ifactual);
        
        %somando 1 LEF
        ifacelef1=parameter(1,3,ifactual);
        termo0=mulef*norma*parameter(1,1,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef1,nflag,lef,rel,parameter,weightDMP,termo0);
        %somando 2 LEF
        ifacelef2=parameter(1,4,ifactual);
        termo0=mulef*norma*parameter(1,2,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacelef2,nflag,lef,rel,parameter,weightDMP,termo0);
        
        % somando 1 REL
        ifacerel1= parameter(2,3,ifactual);
        termo0=murel*norma*parameter(2,1,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel1,nflag,rel,lef,parameter,weightDMP,termo0);
        % somando 2 REL
        ifacerel2= parameter(2,4,ifactual);
        termo0=murel*norma*parameter(2,2,ifactual)*mobility(ifactual);
        [M,I]=auxmatrixassembleDMPopostface(M,I,ifacerel2,nflag,rel,lef,parameter,weightDMP,termo0);
        
    end
    
    x=1000;
    y=1000;
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