function [M,I]=assemblematrixlfvLPEW(parameter,fonte,w,s,nflag,weightDMP,wells,mobility)
global inedge coord bedge bcflag elem elemarea esurn1 esurn2
% incialização das matrizes
I=sparse(size(elem,1),1);
M=sparse(size(elem,1),size(elem,1));
% fonte
I=I+fonte;
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
            % implementação da matriz global no contorno
            M(lef,lef)=M(lef,lef)+ normcont*(parameter(1,1,ifacont)+parameter(1,2,ifacont));
            nolef1=parameter(1,3,ifacont);
            nolef2=parameter(1,4,ifacont);
            % contribuições do nó 1
            if nflag(nolef1,1)<200
                I(lef,1)=I(lef,1)+ normcont*(parameter(1,1,ifacont)*nflag(nolef1,2));
            else
               
                for j=1:(esurn2(nolef1+1)-esurn2(nolef1))
                    
                    post_cont=esurn2(nolef1)+j;
                    auxesurn1=esurn1(post_cont);
                    M(lef, auxesurn1)=M(lef, auxesurn1)-normcont*parameter(1,1,ifacont)*w(post_cont);
                end
            end
            % contribuições do nó 2
            if nflag(nolef2,1)<200
                
                I(lef,1)=I(lef,1)+ normcont*parameter(1,2,ifacont)*nflag(nolef2,2);
                
            else
                
                for j=1:(esurn2(nolef2+1)-esurn2(nolef2))
                    
                    post_cont=esurn2(nolef2)+j;
                    auxesurn1=esurn1(post_cont);
                    M(lef, auxesurn1)=M(lef, auxesurn1)-normcont*parameter(1,2,ifacont)*w(post_cont);
                end
            end
            
        end
    end
end
% Montagem da matriz global

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
    % os nós que conforman os pontos de interpolação no elemento a esquerda
    auxnolef1=parameter(1,3,ifactual);
    auxnolef2=parameter(1,4,ifactual);
    % os nós que conforman os pontos de interpolação no elemento a direita
    auxnorel1=parameter(2,3,ifactual);
    auxnorel2=parameter(2,4,ifactual);
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*murel*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    ARR=norma*mulef*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ ALL;
    M(lef,rel)=M(lef,rel)- ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ ARR;
    M(rel,lef)=M(rel,lef)- ALL;
    % contribuições esquerda
    % termo 1
    if nflag(auxnolef1,1)>200
        
        for j=1:(esurn2(auxnolef1+1)-esurn2(auxnolef1))
            
            post_cont=esurn2(auxnolef1)+j;
            auxesurn1=esurn1(post_cont);
            M(lef, auxesurn1)=M(lef, auxesurn1)-murel*norma*parameter(1,1,ifactual)*w(post_cont);
            M(rel, auxesurn1)=M(rel, auxesurn1)+murel*norma*parameter(1,1,ifactual)*w(post_cont);
        end
    else
        I(lef,1)=I(lef,1)+ murel*norma*parameter(1,1,ifactual)*nflag(auxnolef1,2);
        I(rel,1)=I(rel,1)- murel*norma*parameter(1,1,ifactual)*nflag(auxnolef1,2);
    end
    % no contorno de Neumann com fluxo diferente de zeros
    if nflag(auxnolef1,1)==202
        
        I(lef)=I(lef)+murel*norma*parameter(1,1,ifactual)*s(auxnolef1); %ok
        
        I(rel)=I(rel)-murel*norma*parameter(1,1,ifactual)*s(auxnolef1); %ok
    end
    % termo 2
    if nflag(auxnolef2,1)>200
        
        for j=1:(esurn2(auxnolef2+1)-esurn2(auxnolef2))
            
            post_cont=esurn2(auxnolef2)+j;
            auxesurn1=esurn1(post_cont);
            M(lef, auxesurn1)=M(lef, auxesurn1)- murel*norma*parameter(1,2,ifactual)*w(post_cont);
            M(rel, auxesurn1)=M(rel, auxesurn1)+ murel*norma*parameter(1,2,ifactual)*w(post_cont);
        end
    else
        I(lef,1)=I(lef,1)+ murel*norma*parameter(1,2,ifactual)*nflag(auxnolef2,2);
        I(rel,1)=I(rel,1)- murel*norma*parameter(1,2,ifactual)*nflag(auxnolef2,2);
        
    end
    % no contorno de Neumann com fluxo diferente de zeros
    if nflag(auxnolef2,1)==202
        
        I(lef)=I(lef)+murel*norma*parameter(1,2,ifactual)*s(auxnolef2); %ok
        
        I(rel)=I(rel)-murel*norma*parameter(1,2,ifactual)*s(auxnolef2); %ok
    end
    % contribuição do elemento a direita
    % termo 1
    if nflag(auxnorel1,1)>200
        
        for j=1:(esurn2(auxnorel1+1)-esurn2(auxnorel1))
            
            post_cont=esurn2(auxnorel1)+j;
            auxesurn1=esurn1(post_cont);
            M(lef, auxesurn1)=M(lef, auxesurn1)+mulef*norma*parameter(2,1,ifactual)*w(post_cont);
            M(rel, auxesurn1)=M(rel, auxesurn1)-mulef*norma*parameter(2,1,ifactual)*w(post_cont);
        end
    else
        I(lef,1)=I(lef,1)- mulef*norma*parameter(2,1,ifactual)*nflag(auxnorel1,2);
        I(rel,1)=I(rel,1)+ mulef*norma*parameter(2,1,ifactual)*nflag(auxnorel1,2);
    end
    % no contorno de Neumann com fluxo diferente de zeros
    if nflag(auxnorel1,1)==202
        
        I(lef)=I(lef)-mulef*norma*parameter(2,1,ifactual)*s(auxnorel1); %ok
        
        I(rel)=I(rel)+mulef*norma*parameter(2,1,ifactual)*s(auxnorel1); %ok
    end
    % termo 2
    if nflag(auxnorel2,1)>200
        
         for j=1:(esurn2(auxnorel2+1)-esurn2(auxnorel2))
             
             post_cont=esurn2(auxnorel2)+j;
             auxesurn1=esurn1(post_cont);
             M(lef, auxesurn1)=M(lef, auxesurn1)+ mulef*norma*parameter(2,2,ifactual)*w(post_cont);
             M(rel, auxesurn1)=M(rel, auxesurn1)- mulef*norma*parameter(2,2,ifactual)*w(post_cont);
         end
    else
        I(lef,1)=I(lef,1)- mulef*norma*parameter(2,2,ifactual)*nflag(auxnorel2,2);
        I(rel,1)=I(rel,1)+ mulef*norma*parameter(2,2,ifactual)*nflag(auxnorel2,2);
    end
    % no contorno de Neumann com fluxo diferente de zeros
    if nflag(auxnorel2,1)==202
        
        I(lef)=I(lef)-mulef*norma*parameter(2,2,ifactual)*s(auxnorel2); %ok
        
        I(rel)=I(rel)+mulef*norma*parameter(2,2,ifactual)*s(auxnorel2); %ok
    end
end
%% malha 23x23
%  M(357,:)=0*M(357,:);
%  M(357,357)=1;
%  I(357)=1;
%  M(173,:)=0*M(173,:);
%  M(173,173)=1;
%  I(173)=0;
%% malha 11x11
% M(83,:)=0*M(83,:);
% M(83,83)=1;
% I(83)=1;
% M(39,:)=0*M(39,:);
% M(39,39)=1;
% I(39)=0;
%% Malha 31x31
%  M(698,:)=0*M(698,:);
%  M(698,698)=1;
%  I(698)=1;
%  M(264,:)=0*M(264,:);
%  M(264,264)=1;
%  I(264)=0;
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

end