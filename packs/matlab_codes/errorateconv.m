function errorateconv(solanal, p, velanal,flowrate)
global bedge inedge coord elem elemarea benchmark erromethod

% calculo das pressões máximas e minimas
disp('>> Maximum pressure')
pressuremax= max(p)
disp('>> Minimum pressure')
pressuremin= min(p)

switch benchmark
    case 'edqueiroz'
        disp('>> nao precissa calcular erro para este problema!!!')
        
    otherwise
        
        % recupera as velocidades numéricas
        for iface=1:size(bedge,1)+size(inedge,1)
            if iface<size(bedge,1) ||iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2); 
            end
            norma(iface,1)=norm(coord(v2,:)-coord(v1,:));
            switch benchmark
                case 'lepotier'
                    velnum(iface,1)=flowrate(iface);
                otherwise
                    velnum(iface,1)=flowrate(iface)/norma(iface,1);
            end
            
        end
        switch erromethod
            case 'erromethod1'
                %% O calculo destes erros foram adaptados de Gao an Wu 2010.
                % calcula o erro respeito a pressão
                s=0;
                for i=1:size(elem,1)
                    s=s+(solanal(i,1)-p(i,1))^2*elemarea(i,1);
                end
                
                erropressure=sqrt(s/sum(elemarea));
                % calcula o erro respeito a velocidade
                Q=zeros(size(inedge,1),1);
                for i=1:size(inedge,1)+size(bedge,1)
                    if i>size(bedge,1)
                        Q(i,1)=elemarea(inedge(i-size(bedge,1),3))+elemarea(inedge(i-size(bedge,1),4));
                    else
                        Q(i,1)=elemarea(bedge(i,3));
                    end
                end
                e=-velanal-velnum;
                er=e.^2;
                errovelocity=sqrt((Q'*er)/sum(Q'));
            case 'erromethod2'
                %% O calculo destes erros foram adaptados de Lipnikov et al 2010
                % calcula o erro respeito a pressão
                s1=0;
                s2=0;
                for i=1:size(elem,1)
                    s1=s1+(solanal(i,1)-p(i,1))^2*elemarea(i,1);
                    s2=s2+solanal(i,1)^2*elemarea(i,1);
                end
                erropressure=sqrt(s1/s2);
                % calcula o erro respeito a velocidade
                Q=zeros(size(inedge,1),1);
                for i=1:size(inedge,1)+size(bedge,1)
                    if i>size(bedge,1)
                        Q(i,1)=0.5*(elemarea(inedge(i-size(bedge,1),3))+elemarea(inedge(i-size(bedge,1),4)));
                    else
                        Q(i,1)=elemarea(bedge(i,3));
                    end
                end
                e=-velanal-velnum;
                errovelocity=sqrt((Q'*e.^2)/sum(velanal.^2.*Q));
            case 'erromethod3'
                %% O calculo destes erros foram adaptados de Eigestad 2005
                % calcula o erro respeito a pressão
                s1=0;
                
                for i=1:size(elem,1)
                    s1=s1+(solanal(i,1)-p(i,1))^2*elemarea(i,1);
                end
                erropressure=sqrt(s1);
                % calcula o erro respeito a velocidade
                Q=zeros(size(inedge,1),1);
                for i=1:size(inedge,1)+size(bedge,1)
                    if i>size(bedge,1)
                        Q(i,1)=0.5*(elemarea(inedge(i-size(bedge,1),3))+elemarea(inedge(i-size(bedge,1),4)));
                    else
                        Q(i,1)=elemarea(bedge(i,3));
                    end
                end
                e=-velanal-velnum;
                errovelocity=sqrt((Q'*e.^2)/sum(velanal.^2.*Q));
            case 'erromethod4'
                %% O calculo destes erros foram adaptados de Sheng e Yuan 2015
                % calcula o erro respeito a pressão
                s1=0;
                
                for i=1:size(elem,1)
                    s1=s1+(solanal(i,1)-p(i,1))^2*elemarea(i,1);
                end
                erropressure=sqrt(s1);
                % calcula o erro respeito a velocidade
                Q=zeros(size(inedge,1),1);
                for i=1:size(inedge,1)+size(bedge,1)
                    if i>size(bedge,1)
                        Q(i,1)=0.5*(elemarea(inedge(i-size(bedge,1),3))+elemarea(inedge(i-size(bedge,1),4)));
                    else
                        Q(i,1)=elemarea(bedge(i,3));
                    end
                end
                e=-velanal-velnum;
                errovelocity=sqrt(sum(e.^2));
            case 'erromethod5'
                s1=0;
                
                for i=1:size(elem,1)
                    s1=s1+(solanal(i,1)-p(i,1))^2*elemarea(i,1);
                end
                erropressure=sqrt(s1);
                % calcula o erro respeito a velocidade
                e=-velanal-velnum;
                errovelocity=sqrt(sum(e.^2));
            case 'erromethod6'
                
                % Equacao (42) do artigo Starnoni
               
                s1=((solanal(:,1)-p(:,1)).^2).*elemarea(:,1);
                s2=((solanal(:,1)).^2).*elemarea(:,1);
              
                erropressure=sqrt(sum(s1))/sqrt(sum(s2));
                % calcula o erro respeito a velocidade
                e=-velanal-velnum;
                % Equacao (43) do artigo Starnoni
                errovelocity=sqrt(sum((norma.^2).*e.^2));  
                
        end
        disp('>> Error level in the pressure field')
        erropressure
        disp('>>  Error level in the velocity field')
        errovelocity
end
end