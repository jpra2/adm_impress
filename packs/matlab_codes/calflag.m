function nflag= calflag(pmetodo)
global bedge bcflag coord

switch pmetodo
    case {'nlfvDMPSY','lfvHP','nlfvDMPV1'}
        % determinar o flag do nó interior e fronteira de Neumann
        
        %nflag=ones(size(bedge,1),2); % flags 7.15934 é uma constante
        
        for ifacont=1:size(bedge,1)
            x=bcflag(:,1)==bedge(ifacont,5);
            r=find(x==1);
            nflag(ifacont,2)=bcflag(r,2);
            nflag(ifacont,1)=bcflag(r,1);
            
        end
    case {'mpfad','lfvLPEW','nlfvLPEW','tpfa'}
        nflag=50000*ones(size(coord,1),2); % flags 
        for ifacont=1:size(bedge,1)
            
            x=bcflag(:,1)==bedge(ifacont,4);
            r=find(x==1);
            nflag(bedge(ifacont,1),2)=bcflag(r,2);
            nflag(bedge(ifacont,1),1)=bcflag(r,1);
            
        end
end
end