function [nflag]= contflagface
global bcflag coord elem benchmark bedge

% determinar o flag do nó interior e fronteira de Neumann

%nflag=ones(size(bedge,1),2); % flags 7.15934 é uma constante

for ifacont=1:size(bedge,1)
    a=(coord(bedge(ifacont,1),:)+ coord(bedge(ifacont,2),:))*0.5;
    % coordenadas dos vertices do contorno
    x=a(1,1);
    y=a(1,2);
    % elemento
    lef=bedge(ifacont,3);
    switch benchmark
        
        case {'miao'}
            xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            
            if x<=0.5
               nflag(ifacont,2)=14*x+y;
            else
              nflag(ifacont,2)= 4*x+y+5;
              
            end
        case {'starnonigrav1'}
             xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            
            if nflag(ifacont,1)>200
               nflag(ifacont,2)=0;
            else
                if single(y)>=0.5
                    nflag(ifacont,2)= 11-10*y;
                else
                    
                    % condicao de contorno de Dirichlet no lado direito e esquerdo
                    nflag(ifacont,2)= 6.5-1*y;
                end
                 
            end
        case{'starnonigrav2'}
            xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            
            if nflag(ifacont,1)>200
               nflag(ifacont,2)=0;
            else
               nflag(ifacont,2)=1+sin(x)*cos(y);  
            end
        case {'starnonigrav3'}
             xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            
            if nflag(ifacont,1)>200
               nflag(ifacont,2)=0;
            else
               if single(y)>0.5
                    nflag(ifacont,2)=sin(x)*cos(y)+ 11-10*y;
                else
                    
                    % condicao de contorno de Dirichlet no lado direito e esquerdo
                    nflag(ifacont,2)=sin(x)*cos(y)+ 6.5-1*y;
                end 
                 
            end 
         case {'starnonigrav4'}
             xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            
            if nflag(ifacont,1)>200
               nflag(ifacont,2)=0;
            else
               if single(y)>0.5
                    nflag(ifacont,2)=100*sin(x)*cos(y)+11- 10*y;
                else
                    
                    % condicao de contorno de Dirichlet no lado direito e esquerdo
                    nflag(ifacont,2)=100*sin(x)*cos(y)+6.5- 1*y;
                end 
                 
            end   
        
        case{'zhangkobaise'}
            xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            teta=calculoteta(x,y);
            alfa=0.5645;
            r=sqrt(x^2+y^2);
            if elem(lef,5)==1
                a1=1.0;
                b1=-12.0414;
                nflag(ifacont,2)=10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
            elseif elem(lef,5)==2
                a2=-4.8591;
                b2=-6.0699;
                nflag(ifacont,2)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
            else
                a3=-0.9664;
                b3=-0.2837;
                nflag(ifacont,2)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
            end
        case {'zhangkobaise2'}
            xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            teta=calculoteta(x,y);
            alfa=0.6142;
            r=sqrt(x^2+y^2);
            if elem(lef,5)==1
                a1=1.0;
                b1=-1.0546;
                nflag(ifacont,2)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
            elseif elem(lef,5)==2
                a2=-0.4275;
                b2=0.2142;
                nflag(ifacont,2)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
            else
                a3=-0.7604;
                b3=-0.6495;
                nflag(ifacont,2)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
            end
        case {'zhangkobaise3'}
            xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            teta=calculoteta(x,y);
            alfa=0.8866;
            r=sqrt(x^2+y^2);
            if elem(lef,5)==1
                a1=1.0;
                b1=-0.3706;
                nflag(ifacont,2)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
            elseif elem(lef,5)==2
                a2=-0.0144;
                b2=0.0022;
                
                nflag(ifacont,2)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
            else
                a3=0.7544;
                b3=-0.6564;
                
                nflag(ifacont,2)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
            end
            
            
        case {'homogeneo', 'heterogeneo','nikitin','durlofsky','lamine','shuec'}
            x=bcflag(:,1)==bedge(ifacont,5);
            r=find(x==1);
            nflag(ifacont,2)=bcflag(r,2);
            nflag(ifacont,1)=bcflag(r,1);
            %%
        case 'crumpton'
            alfa=1000;
            if x<0 || x==0
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=(2*sin(y)+cos(y))*alfa*x + sin(y)+3*alfa;
            else
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=exp(x)*sin(y)+3*alfa;
            end
            
        case 'crumptonhyman'
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=exp(x*y);
        case 'gaowu2'
            %%
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=0.5*((sin((1-x)*(1-y))/(sin(1)))+(1-x)^3*(1-y)^2);
        case 'gaowu1'
            %%
            if (x<0.5 || x==0.5)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=(1+(x-0.5)*(0.1+8*pi*(y-0.5)))*exp(-20*pi*((y-0.5)^2));
                
            else
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=exp(x-0.5)*exp(-20*pi*((y-0.5)^2));
                
            end
            %%
        case  'gaowu3'
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=exp(-20*pi*((x-0.5)^2 + (y-0.5)^2));
        case 'lepotier'
            %%
            xx=bcflag(:,1)==bedge(ifacont,5);
            rr=find(xx==1);
            nflag(ifacont,1)=bcflag(rr,1);
            % no artigo de lipnikov
            %nflag(ifacont,2)=sin(pi*x)*sin(pi*y);
            %no artigo Zhang e Kobaise
            nflag(ifacont,2)=sin(pi*x)*sin(pi*y)+1;
        case 'lipnikov1'
            %%
            if x<0.5 || x==0.5
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=1-2*y^2+4*x*y+6*x+2*y;
                
            else
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=-2*y^2+1.6*x*y-0.6*x+3.2*y+4.3;
                
            end
        case 'edwards'
            %%
            f=(4/((50-2)*0.1+1));
            b2=(0.1-1)*f;
            c1=50*0.1*f;
            c2=f;
            d2=-f*(1/10);
            d1=d2;
            if x<0.5
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=c1*x^2+d1*y^2;
                
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
            if (y==0.4444)||(y==0.5555)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=2;
            elseif (y==0)||(y==1)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=0;
            end
            if (x==0.4444)||(x==0.5555)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=2;
            elseif (x==0)||(x==1)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=0;
            end
        case 'guangwei'
            %%
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=0;
        case 'guangwei1'
            %%
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=0;
        case 'gaowu4'
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=sin(pi*x)*sin(pi*y);
        case 'gaowu5'
            %%
            if (((0<x || x==0 )&& (x <0.2 || x==0.2)) && y==0) || (((0<y || y==0 )&& (y <0.2 || y==0.2)) && x==0)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=1;
            elseif (((0.8<x || x==0.8 )&& (x <1 || x==1)) && y==1) || (((0.8<y || y==0.8 )&& (y <1 || y==1)) && x==1)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=0;
            elseif (((0.3<x || x==0.3 )&& (x <1 || x==1)) && y==0) || (((0.3<y || y==0.3 )&& (y <1 || y==1)) && x==0)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=0.5;
            elseif (((0<x || x==0 )&& (x <0.7 || x==0.7)) && y==1) || (((0<y || y==0 )&& (y <0.7 || y==0.7)) && x==1)
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=0.5;
            else
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=0.5;
            end
        case 'gaowu6'
            %%
            nflag(ifacont,1)=201;
            nflag(ifacont,2)=0;
        case 'gaowu7'
            %%
            delta=0.2;
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=-x-delta*y;
        case 'gaowu8'
            %%
            delta=0.2;
            phi1=y-delta*(x-0.5)-0.475;
            phi2=phi1-0.05;
            % dominio 1
            if phi1<0
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=-phi1;
                %dominio
            elseif phi1>0 && phi2<0
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=-100*phi1;
            elseif phi2>0
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=-phi2-5;
            end
        case 'gaowu9'
            %%
            delta=0.2;
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=2-x-delta*y;
        case {'benchmar5_7','benchmar5_6'}
            %%
            
            nflag(ifacont,1)=101;
            nflag(ifacont,2)=0;
        case 'edqueiroz'
            %%
            if bedge(ifacont,5)==101
                nflag(ifacont,1)=101;
                nflag(ifacont,2)=0;
            else
                nflag(ifacont,1)=102;
                nflag(ifacont,2)=2;
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