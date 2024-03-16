function [M,I] = globalmatrixtpfa(Kde, Kn, nflagface, Hesq,gravresult,gravrate,gravno,gravelem,fonte)

global inedge bedge elem coord bcflag gravitational strategy

% Constrói a matriz global.
% prealocação da matriz global e do vetor termo de fonte
M=zeros(size(elem,1),size(elem,1));
I=zeros(size(elem,1),1);
I=I+fonte;
% Loop de faces de contorno
m=0;

for ifacont=1:size(bedge,1)
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    normcont=norm(v0);
    lef=bedge(ifacont,3);
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
    
    if bedge(ifacont,5)<200
        
        c1=nflagface(ifacont,2);
        
        if strcmp(gravitational,'yes')
            if strcmp(strategy,'starnoni')
                m=gravrate(ifacont);
            else
                % proposto de nos
                m1=-nflagface(ifacont,2);
                m=A*(norm(v0)^2*m1-norm(v0)^2*gravelem(lef));
            end
        else
            m=0;
        end
        %Preenchimento
        M(bedge(ifacont,3),bedge(ifacont,3))=M(bedge(ifacont,3),bedge(ifacont,3))- A*(norm(v0)^2);
        
        I(bedge(ifacont,3))=I(bedge(ifacont,3))-c1*A*(norm(v0)^2)+m;
        
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(bedge(ifacont,3))=I(bedge(ifacont,3))- normcont*bcflag(r,2);
        
        
    end
    
end

for iface=1:size(inedge,1),
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Contabiliza as contribuições do fluxo numa faces  para os elementos %
    %a direita e a esquerda dela.                                        %
    M(inedge(iface,3), inedge(iface,3))=M(inedge(iface,3), inedge(iface,3))-Kde(iface,1);
    M(inedge(iface,3), inedge(iface,4))=M(inedge(iface,3), inedge(iface,4))+Kde(iface,1);
    M(inedge(iface,4), inedge(iface,4))=M(inedge(iface,4), inedge(iface,4))-Kde(iface,1);
    M(inedge(iface,4), inedge(iface,3))=M(inedge(iface,4), inedge(iface,3))+Kde(iface,1);
    if strcmp(gravitational,'yes')
        if strcmp(strategy,'starnoni')
            m=gravrate(size(bedge,1)+iface);
        else
            m= Kde(iface)*(gravelem(rel,1)-gravelem(lef,1));
            %m=gravrate(size(bedge,1)+iface,1);
            
        end
        I(lef)=I(lef)+m;
        I(rel)=I(rel)-m;
    end
        
end
end
