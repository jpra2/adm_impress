% Esta funcao calcula ou impoe todos flags de contorno e seus respectivos
% valores
function nflag= contflagno
global  bcflag coord elem benchmark bedge

% quarta coluna corresponde a flag do vertice
% quinta coluna corresponde a flag da face
% 500000: simplesmente representa o flag dos vertices interiores
nflag=50000*ones(size(coord,1),2); 

for ifacont=1:size(bedge,1)
    a=coord(bedge(ifacont,1),:);
    x=a(1,1);
    y=a(1,2);
    
    lef=bedge(ifacont,3);
    xx=bcflag(:,1)==bedge(ifacont,4); % quarta coluna bedge 
    %corresoponde ao flag do vertice
    rr=find(xx==1);
    nflag(bedge(ifacont,1),1)=bcflag(rr,1);
    switch benchmark
        case 'miao'
  
            if x<=0.5
               nflag(bedge(ifacont,1),2)=14*x+y;
            else
              nflag(bedge(ifacont,1),2)= 4*x+y+5;
              
            end
        case 'starnonigrav1'
            % flag do no
            xx=bcflag(:,1)==bedge(ifacont,4); % quarta coluna bedge 
            %corresoponde ao flag do vertice
            rr=find(xx==1);
            nflag(bedge(ifacont,1),1)=bcflag(rr,1);
             h1=10; h2=1;
            if nflag(bedge(ifacont,1),1)>200
              % condicao de contorno de Neumann no lado superior e inferior
               nflag(bedge(ifacont,1),2)=0;
            else
                % condicao de contorno de Dirichlet no lado direito e esquerdo
                %nflag(bedge(ifacont,1),2)=(11-h1*y)*logical(y>=0.5)+(6.5-h2*y)*logical(y<0.5);
                if y>=0.5
                    nflag(bedge(ifacont,1),2)= 11-h1*y;
                else    
                    nflag(bedge(ifacont,1),2)= 6.5-h2*y;
                end
              
            end
        case 'starnonigrav2'
                        
            if nflag(bedge(ifacont,1),1)>200
              % condicao de contorno de Neumann no lado superior e inferior
               nflag(bedge(ifacont,1),2)=0;
            else
              % condicao de contorno de Dirichlet no lado direito e esquerdo
              
              nflag(bedge(ifacont,1),2)= 1+sin(x)*cos(y);
            end
        case 'starnonigrav3'
            xx=bcflag(:,1)==bedge(ifacont,4);
            rr=find(xx==1);
            % sin e cos calcula o angulo em radianes
            nflag(bedge(ifacont,1),1)=bcflag(rr,1);
            h1=10; h2=1;
            if nflag(bedge(ifacont,1),1)>200
              % condicao de contorno de Neumann no lado superior e inferior
               nflag(bedge(ifacont,1),2)=0;
            else
                if y>=0.5
                    nflag(bedge(ifacont,1),2)= sin(x)*cos(y)+11-h1*y;
                else
                    % condicao de contorno de Dirichlet no lado direito e esquerdo
                    nflag(bedge(ifacont,1),2)= sin(x)*cos(y)+6.5-h2*y;
                end
              
            end
         case 'starnonigrav4'
            xx=bcflag(:,1)==bedge(ifacont,4);
            rr=find(xx==1);
            nflag(bedge(ifacont,1),1)=bcflag(rr,1);
            h1=10; h2=1;
            if nflag(bedge(ifacont,1),1)>200
              % condicao de contorno de Neumann no lado superior e inferior
               nflag(bedge(ifacont,1),2)=0;
            else
                if y>=0.5
                    nflag(bedge(ifacont,1),2)= 100*sin(x)*cos(y)+11-h1*y;
                else
                    
                    % condicao de contorno de Dirichlet no lado direito e esquerdo
                    nflag(bedge(ifacont,1),2)= 100*sin(x)*cos(y)+6.5-h2*y;
                end
            end
       
        case{'zhangkobaise'}
            xx=bcflag(:,1)==bedge(ifacont,4);% quarta coluna bedge 
            %corresoponde ao flag do vertice
            rr=find(xx==1);
            nflag(bedge(ifacont,1),1)=bcflag(rr,1);
            teta=calculoteta(x,y);
            alfa=0.5645;
            r=sqrt(x^2+y^2);
            if elem(lef,5)==1
                a1=1.0;
                b1=-12.0414;
                nflag(bedge(ifacont,1),2)=10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
            elseif elem(lef,5)==2
                a2=-4.8591;
                b2=-6.0699;
                nflag(bedge(ifacont,1),2)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
            else
                a3=-0.9664;
                b3=-0.2837;
                nflag(bedge(ifacont,1),2)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
            end
        case {'zhangkobaise2'}
            xx=bcflag(:,1)==bedge(ifacont,4);
            rr=find(xx==1);
            nflag(bedge(ifacont,1),1)=bcflag(rr,1);
            teta=calculoteta(x,y);
            alfa=0.6142;
            r=sqrt(x^2+y^2);
            if elem(lef,5)==1
                a1=1.0;
                b1=-1.0546;
                nflag(bedge(ifacont,1),2)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
            elseif elem(lef,5)==2
                a2=-0.4275;
                b2=0.2142;
                nflag(bedge(ifacont,1),2)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
            else
                a3=-0.7604;
                b3=-0.6495;
                nflag(bedge(ifacont,1),2)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
            end
        case {'zhangkobaise3'}
            xx=bcflag(:,1)==bedge(ifacont,4);
            rr=find(xx==1);
            nflag(bedge(ifacont,1),1)=bcflag(rr,1);
            teta=calculoteta(x,y);
            alfa=0.8866;
            r=sqrt(x^2+y^2);
            if elem(lef,5)==1
                a1=1.0;
                b1=-0.3706;
                nflag(bedge(ifacont,1),2)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
            elseif elem(lef,5)==2
                a2=-0.0144;
                b2=0.0022;
                
                nflag(bedge(ifacont,1),2)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
            else
                a3=0.7544;
                b3=-0.6564;
                
                nflag(bedge(ifacont,1),2)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
            end
            
        case {'homogeneo', 'heterogeneo','nikitin','durlofsky','lamine',...
                'shuec','altamenteheterogeneo','nikitin1','pinchout'}
            x=bcflag(:,1)==bedge(ifacont,4);
            r=find(x==1);
            nflag(bedge(ifacont,1),2)=bcflag(r,2);
            nflag(bedge(ifacont,1),1)=bcflag(r,1);
        case 'zigzagfract'
            %%
            x=bcflag(:,1)==bedge(ifacont,4);
            r=find(x==1);
            nflag(bedge(ifacont,1),2)=bcflag(r,2);
            nflag(bedge(ifacont,1),1)=bcflag(r,1);
        case 'crumpton'
            
            %%
            alfa=1000;
            if x<0 || x==0
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=(2*sin(y)+cos(y))*alfa*x + sin(y)+3*alfa;
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=exp(x)*sin(y)+3*alfa;
            end
        case 'crumptonhyman'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=exp(x*y);
        case 'gaowu2'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0.5*((sin((1-x)*(1-y))/(sin(1)))+(1-x)^3*(1-y)^2);
            
        case 'gaowu1'
            %%
            if (x<0.5 || x==0.5)
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=(1+(x-0.5)*(0.1+8*pi*(y-0.5)))*exp(-20*pi*((y-0.5)^2));
                
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=exp(x-0.5)*exp(-20*pi*((y-0.5)^2));
                
            end
            
        case  'gaowu3'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=exp(-20*pi*((x-0.5)^2 + (y-0.5)^2));
        case 'gaowu4'
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=sin(pi*x)*sin(pi*y);
        case 'lepotier'
            %%
            
            nflag(bedge(ifacont,1),1)=101;
            % no artigo de Lipnikov
            %nflag(bedge(ifacont,1),2)=sin(pi*x)*sin(pi*y);
            % no artigo de Zhang e Kobaise
            nflag(bedge(ifacont,1),2)=sin(pi*x)*sin(pi*y);
        case 'lipnikov1'
            %%
            if x<0.5 || x==0.5
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=1-2*y^2+4*x*y+6*x+2*y;
                
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=-2*y^2+1.6*x*y-0.6*x+3.2*y+4.3;
                
            end
        case 'edwards'
            %%
            a1=1;
            f=((4*a1)/(((1/50)-2)*(1/10)+1));
            b2=((1/10)-1)*f;
            c2=f;
            d2=-c2*(1/10);
            c1=(1/50)*(1/10)*f;
            d1=d2;
            if x<0.5
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=c1*x^2+d1*y^2;
                
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=1+b2*x+c2*x^2+d2*y^2;
                
            end
        case 'shenyuan16'
            
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=sin(pi*x)*sin(pi*y);
        case 'herbin'
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=sin(pi*x)*sin(pi*y);
        case 'herbinhubert'
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=16*x*(1-x)*y*(1-y);
        case 'lipnikov2'
            %%
            %             if (y==0.4444)||(y==0.5555)
            %                 nflag(bedge(ifacont,1),1)=102;
            %                 nflag(bedge(ifacont,1),2)=12;
            %             elseif (y==0)||(y==1)
            %                 nflag(bedge(ifacont,1),1)=101;
            %                 nflag(bedge(ifacont,1),2)=10;
            %             end
            %             if (x==0.4444)||(x==0.5555)
            %                 nflag(bedge(ifacont,1),1)=102;
            %                 nflag(bedge(ifacont,1),2)=12;
            %             elseif (x==0)||(x==1)
            %                 nflag(bedge(ifacont,1),1)=101;
            %                 nflag(bedge(ifacont,1),2)=10;
            %             end
            nflag(bedge(ifacont,1),1)=201;
            nflag(bedge(ifacont,1),2)=0;
        case 'FriisEdwards'
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0;
        case 'lipnikov3'
            %%
            if bedge(ifacont,5)==102
                nflag(bedge(ifacont,1),1)=102;
                nflag(bedge(ifacont,1),2)=2;
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0;
            end
            
        case 'guangwei'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0;
        case 'guangwei1'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0;
            
        case 'gaowu5'
            %%
             if (((0<x || x==0 )&& (x <0.2 || x==0.2)) && y==0) || (((0<y || y==0 )&& (y <0.2 || y==0.2)) && x==0)
                 nflag(bedge(ifacont,1),1)=101;
                 nflag(bedge(ifacont,1),2)=1;
             elseif (((0.8<x || x==0.8 )&& (x <1 || x==1)) && y==1) || (((0.8<y || y==0.8 )&& (y <1 || y==1)) && x==1)
                 nflag(bedge(ifacont,1),1)=101;
                 nflag(bedge(ifacont,1),2)=0;
             elseif (((0.3<x || x==0.3 )&& (x <1 || x==1)) && y==0) || (((0.3<y || y==0.3 )&& (y <1 || y==1)) && x==0)
                 nflag(bedge(ifacont,1),1)=101;
                 nflag(bedge(ifacont,1),2)=0.5;
             elseif (((0<x || x==0 )&& (x <0.7 || x==0.7)) && y==1) || (((0<y || y==0 )&& (y <0.7 || y==0.7)) && x==1)
                 nflag(bedge(ifacont,1),1)=101;
                 nflag(bedge(ifacont,1),2)=0.5;
             else
                 nflag(bedge(ifacont,1),1)=101;
                 nflag(bedge(ifacont,1),2)=0.5;
             end
             
%             nflag(bedge(ifacont,1),1)=101;
%             nflag(bedge(ifacont,1),2)=50;
        case 'gaowu6'
            %%
            nflag(bedge(ifacont,1),1)=201;
            nflag(bedge(ifacont,1),2)=0;
            
        case 'gaowu7'
            %%
            delta=0.2;
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=-x-delta*y;
        case 'gaowu8'
            %%
            delta=0.2;
            phi1=y-delta*(x-0.5)-0.475;
            phi2=phi1-0.05;
            % dominio 1
            if phi1<0 || phi1==0
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=-phi1;
                %dominio
            elseif phi1>0 && phi2<0
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=-100*phi1;
            elseif phi2>0 || phi2==0
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=-phi2-5;
            end
            
            
        case 'gaowu9'
            %%
            delta=0.2;
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=2-x-delta*y;
            
        case {'benchmar5_7','benchmar5_6'}
            %%
            
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0;
            
        case 'edqueiroz'
            %%
            if bedge(ifacont,4)==101
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0;
            else
                nflag(bedge(ifacont,1),1)=102;
                nflag(bedge(ifacont,1),2)=2;
            end 
    end
end
end

function teta=calculoteta(x,y)

if y<0
    if x>0
        m=atan(y/x);
    elseif x<0 && y>=0
        m=atan(y/x)+pi;
    elseif x<0 && y<0
        m=atan(y/x)-pi;
    elseif x==0 && y>0
        m=pi/2;
    elseif x==0 && y<0
        m=-pi/2;
    else
        m='indefined';
    end
    teta=2*pi+m;
else
    if x>0
        m=atan(y/x);
    elseif x<0 && y>=0
        m=atan(y/x)+pi;
    elseif x<0 && y<0
        m=atan(y/x)-pi;
    elseif x==0 && y>0
        m=pi/2;
    elseif x==0 && y<0
        m=-pi/2;
    else
        m='indefined';
    end
    teta=m;
    
end

end