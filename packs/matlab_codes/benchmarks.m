% este codigo contem todas as informacoes dos diferentes casos, como por
% exemplo: tensor de permeabilidade (kmap), pressão analitica (u), termos de fonte (fonte),
% velocidade analitica (vel),gravidade (grav)
function[elem,kmap,normKmap,u,bedge,fonte,vel,grav,gravno,gravface,grav_elem_escalar]=benchmarks(kmap,elem,bedge)
global centelem coord inedge normals elemarea bcflag benchmark
normKmap=0;
vel=0;
u=0;
fonte=0;
grav=zeros(size(elem,1),3);
gravno=0;
gravface=0;
grav_elem_escalar=0;

switch benchmark
    
    case 'miao'
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            if x<=0.5
                u(i,1)= 14*x+y;
                kmap(i,:) = [i 3 1 1 3];
                
            else
                u(i,1)= 4*x+y+5;
                kmap(i,:) = [i 10 3 3 10];
            end
            elem(i,5)=i;
        end
        R1=[0 -1 0; 1 0 0; 0 0 0];
        for iface=1:size(bedge,1)+size(inedge,1)
            
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
            end
            IJ=coord(v2,:)-coord(v1,:);
            norma=norm(IJ);
            nij=R1*IJ'/norma;
            p1=(coord(v2,:)+coord(v1,:))*0.5;
            if p1(1,1)<=0.5
                a=[-43, -17, 0];
                
            else
                a=[-43, -22, 0];
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'starnonigrav1'
        
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            y = centelem(i,2);
            % parametro segundo  Starnoni
            h1=10; h2=1;
            if y>=0.5
                % solucao analitica
                u(i,1)= 11-h1*y;
                
                % calculo do gravidade
                grav(i,:)=h1*[0,1,0];
                grav_elem_escalar(i,1)=(-11+h1*y);
            else
                % solucao analitica
                u(i,1)= 6.5-h2*y;
                % calculo do gravidade
                grav(i,:)=h2*[0,1,0];
                grav_elem_escalar(i,1)=(-6.5+h2*y);
            end
        end
        for jj=1:size(coord,1)
            %Define "x" and "y"
            h1=10; h2=1;
            y2 = coord(jj,2);
            % parametro segundo  Starnoni
            
            if y2>=0.5
                % solucao analitica
                gravno(jj,1)= -11+h1*y2;
            else
                % solucao analitica
                gravno(jj,1)= -6.5+h2*y2;
            end
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                
                ym=a(1,2);
               
                    if ym>=0.5
                        h1=10;
                        V=-[-0.1*0 -0 0];
                    else
                        h2=1;
                        V=-[-0.1*0 -0 0];
                    end
                
            %Obtain the flow rate
                
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                
                ym=a(1,2);
                if ym>=0.5
                    h1=10;
                V=-[-0.1*0 0 0];
                else
                    h2=1;
                V=-[-0.1*0 -0 0];    
                end
            end
            F(j,1) = dot(V,nij');
            y11=a(1,2);
            % parametro segundo  Starnoni
            
            if y11>=0.5
                % solucao analitica
                gravface(j,1:3)= h1*[0,1,0];
            else
                % solucao analitica
                gravface(j,1:3)= h2*[0,1,0];
            end
        end

        vel=F;
        K=kmap;
        elem(:,5)=1;
    case 'starnonigrav2'
         R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            % solucao analitica foi calculado usando pag. 385
            % Calculo II Tom Apostol
            %  sin e cos no radianes
            u(i,1)=1+sin(x)*cos(y);
            
            % gravidade
            grav(i,:)=[-cos(x)*cos(y) sin(x)*sin(y)];
           grav_elem_escalar(i,1)=-1-sin(x)*cos(y);
        end
        for j=1:size(coord,1)
            %Define "x" and "y"
            x1 = coord(j,1);
            y1 = coord(j,2);
            % parametro segundo  Starnoni
            
            % solucao analitica
            gravno(j,1)=-1- sin(x1)*cos(y1);
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
            end
            a=0.5*(coord(v1,:)+coord(v2,:));
            x2=a(1,1);
            y2=a(1,2);
            gravface(j,1:2)= [-cos(x2)*cos(y2) sin(x2)*sin(y2) ];
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);  
            end
            norma=norm(IJ);
            nij=R*IJ'/norma;
            % note que a velocidade analitica: velocidade_pressao +
            % velocidade_gravitacional, embora essa soma foi zero porque
            % g=-grad(p).
            V=-[cos(a(1,1))*cos(a(1,2))-0.1*sin(a(1,1))*sin(a(1,2)),...
                   0.1*cos(a(1,1))*cos(a(1,2))-sin(a(1,1))*sin(a(1,2)),0]+...
                   [cos(a(1,1))*cos(a(1,2))-0.1*sin(a(1,1))*sin(a(1,2)),...
                   0.1*cos(a(1,1))*cos(a(1,2))-sin(a(1,1))*sin(a(1,2)),0];
            F(j,1) = dot(V,nij'); 
        end

        vel=F;
        K=kmap;
        elem(:,5)=1;
    case 'starnonigrav3'
        R=[0 1 0; -1 0 0; 0 0 0];
        h1=10;
        h2=1;
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            % parametro segundo  Starnoni
            
            if single(y)>0.5
                
                % solucao analitica
                u(i,1)= sin(x)*cos(y)+11-h1*y;
                
                % calculo do gravidade
                grav(i,:)=[-cos(x)*cos(y) h1+sin(x)*sin(y)];
            else
                % solucao analitica
                u(i,1)= sin(x)*cos(y)+6.5-h2*y;
                % calculo do gravidade
                grav(i,:)=[-cos(x)*cos(y) h2+sin(x)*sin(y)];
                
            end
        end
        for j=1:size(coord,1)
            %Define "x" and "y"
            x21=coord(j,1);
            y21 = coord(j,2);
            
            if single(y21)>0.5
                
                % solucao analitica
                gravno(j,1)= sin(x21)*cos(y21)-h1*y21;
            else
                % solucao analitica
                gravno(j,1)= sin(x21)*cos(y21)-h2*y21;
                
            end
            
        end
        
        for jj=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if jj<=size(bedge,1)
                v1=bedge(jj,1);
                v2=bedge(jj,2);
            else
                v1=inedge(jj-size(bedge,1),1);
                v2=inedge(jj-size(bedge,1),2);
            end
            a=0.5*(coord(v1,:)+coord(v2,:));
            x11=a(1,1);
            y11=a(1,2);
           
            if single(y11)>0.5
                
                % solucao analitica
                aaa= [-cos(x11)*cos(y11), sin(x11)*sin(y11)+h1];
            else
                % solucao analitica
               aaa= [-cos(x11)*cos(y11), sin(x11)*sin(y11)+h2];
                
            end
            gravface(jj,1:2)=aaa;
        end
         for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);  
            end
            y111=a(1,2);
            if single(y111)>0.5
                V=-[cos(a(1,1))*cos(a(1,2))-0.1*(sin(a(1,1))*sin(a(1,2))+h1),...
                   0.1*cos(a(1,1))*cos(a(1,2))-(sin(a(1,1))*sin(a(1,2))+h1),0]+...
                   [cos(a(1,1))*cos(a(1,2))-0.1*(sin(a(1,1))*sin(a(1,2))+h1),...
                   0.1*cos(a(1,1))*cos(a(1,2))-(sin(a(1,1))*sin(a(1,2))+h1),0];
            else
               V=-[cos(a(1,1))*cos(a(1,2))-0.1*(sin(a(1,1))*sin(a(1,2))+h2),...
                   0.1*cos(a(1,1))*cos(a(1,2))-(sin(a(1,1))*sin(a(1,2))+h2),0]+...
                   [cos(a(1,1))*cos(a(1,2))-0.1*(sin(a(1,1))*sin(a(1,2))+h2),...
                   0.1*cos(a(1,1))*cos(a(1,2))-(sin(a(1,1))*sin(a(1,2))+h2),0]; 
            end
            norma=norm(IJ);
            nij=R*IJ'/norma;
            % note que a velocidade analitica: velocidade_pressao +
            % velocidade_gravitacional, embora essa soma sera zero porque
            % g=-grad(p).
            
            F(j,1) = dot(V,nij'); 
        end

        vel=F;
        K=kmap;
        elem(:,5)=1;
    case 'starnonigrav4'
        % parametro segundo  Starnoni
            h1=10;
            h2=1;
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            
            if single(y)>0.5
                
                % solucao analitica
                u(i,1)= 100*sin(x)*cos(y)+11-h1*y;
               
                % calculo do gravidade
                grav(i,:)=[-100*cos(x)*cos(y) h1+100*sin(x)*sin(y)];
            else
                % solucao analitica
                u(i,1)= 100*sin(x)*cos(y)+6.5-h2*y;
                % calculo do gravidade
                grav(i,:)=[-100*cos(x)*cos(y) h2+100*sin(x)*sin(y)];
                
            end
        end
        for j=1:size(coord,1)
            %Define "x" and "y"
            x21=coord(j,1);
            y21 = coord(j,2);
            
            if single(y21)>0.5
                
                % solucao analitica
                gravno(j,1)= 100*sind(x21)*cosd(y21)-h1*y21;
            else
                % solucao analitica
                gravno(j,1)= 100*sind(x21)*cosd(y21)-h2*y21;
                
            end
            
        end
        
        for jj=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if jj<=size(bedge,1)
                v1=bedge(jj,1);
                v2=bedge(jj,2);                
            else
                v1=inedge(jj-size(bedge,1),1);
                v2=inedge(jj-size(bedge,1),2);
            end
            a=0.5*(coord(v1,:)+coord(v2,:));
            x11=a(1,1);
            y11=a(1,2);
            
            if single(y11)>0.5
               bbb= [-100*cos(x11)*cos(y11), 100*sin(x11)*sin(y11)+h1];
            else
               bbb= [-100*cos(x11)*cos(y11), 100*sin(x11)*sin(y11)+h2];
                
            end           
           gravface(jj,1:2)=bbb; 
        end
        
        K=kmap;
        elem(:,5)=1;
    case 'zhangkobaise'
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            teta=calculoteta(x,y);
            alfa=0.5645;
            r=sqrt(x^2+y^2);
            coef=alfa*(alfa-1)*(x^2+y^2)^(0.5*alfa-2);
            
            if elem(i,5)==1
                a1=1.0;
                b1=-12.0414;
                
                u(i,1)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
                kmap(1,:) = [1 1 0 0 1];
                fonte(i,1) =0;
            elseif elem(i,5)==2
                a2=-4.8591;
                b2=-6.0699;
                
                u(i,1)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                kmap(2,:) = [2 10 0 0 10];
                %==========================================================
                sum1= (x^2)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                sum2=(x*y)*(-2*b2*cos(alfa*teta)+2*a2*sin(alfa*teta));
                sum3=(y^2)*(-a2*cos(alfa*teta)-b2*sin(alfa*teta));
                sumtotal=sum1+sum2+sum3;
                fonte(i,1) =(10-10)*coef*sumtotal*elemarea(i,1);
            else
                a3=-0.9664;
                b3=-0.2837;
                
                u(i,1)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
                kmap(3,:) = [3 100 0 0 100];
                fonte(i,1) =0;
            end
            
            
        end  %End of FOR
        K=kmap;
        
        %% velocities
        for ianalit = 1:size(bedge,1)
            lef=bedge(ianalit,3);
            IJ=coord(bedge(ianalit,2),:)-coord(bedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            
            auxpoint=(coord(bedge(ianalit,2),:)+coord(bedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            teta1=calculoteta(x,y);
            alfa=0.5645;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-12.0414;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-4.8591;
                b2=-6.0699;
                k11=10;
                k22=10;
                % dtheta/dx
                dtetadx=calculotetadx(x,y);
                %==========================================================
                var1=((alfa*(r^alfa))/r^2)*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                var2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1));
                %
                part1=var1*x;
                part2=var2*dtetadx;
                % dtheta/dy
                dtetady=calculotetady(x,y);
                part3=var1*y;
                part4=var2*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=-0.9664;
                b3=-0.2837;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            %Obtain the flow rate
            F(ianalit,1) = dot(V,nij');
        end  %End of FOR (boundary edges)
        
        for ianalit=1:size(inedge,1)
            lef=inedge(ianalit,3);
            
            IJ=coord(inedge(ianalit,2),:)-coord(inedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(ianalit,2),:)+coord(inedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            
            teta1=calculoteta(x,y);
            alfa=0.5645;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-12.0414;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-4.8591;
                b2=-6.0699;
                k11=10;
                k22=10;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=-0.9664;
                b3=-0.2837;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            
            F(ianalit+size(bedge,1),1) = dot(V,nij');
        end
        
        vel=F;
    case 'zhangkobaise2'
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            teta=calculoteta(x,y);
            alfa=0.6142;
            r=sqrt(x^2+y^2);
            coef=alfa*(alfa-1)*(x^2+y^2)^(0.5*alfa-2);
            if elem(i,5)==1
                a1=1.0;
                b1=-1.0546;
                
                u(i,1)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
                kmap(1,:) = [1 1 0 0 1];
                fonte(i,1) =0;
            elseif elem(i,5)==2
                a2=-0.4275;
                b2=0.2142;
                
                u(i,1)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                kmap(2,:) = [2 10 0 0 1000];
                %==========================================================
                sum1= (x^2)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                sum2=(x*y)*(-2*b2*cos(alfa*teta)+2*a2*sin(alfa*teta));
                sum3=(y^2)*(-a2*cos(alfa*teta)-b2*sin(alfa*teta));
                sumtotal=sum1+sum2+sum3;
                fonte(i,1) =(10-1000)*coef*sumtotal*elemarea(i,1);
            else
                a3=-0.7604;
                b3=-0.6495;
                
                u(i,1)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
                kmap(3,:) = [3 100 0 0 100];
                fonte(i,1) =0;
            end
            
            
        end  %End of FOR
        K=kmap;
        
        %% velocities
        for ianalit = 1:size(bedge,1)
            lef=bedge(ianalit,3);
            IJ=coord(bedge(ianalit,2),:)-coord(bedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            
            auxpoint=(coord(bedge(ianalit,2),:)+coord(bedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            teta1=calculoteta(x,y);
            alfa=0.6142;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-1.0546;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-0.4275;
                b2=0.2142;
                k11=10;
                k22=1000;
                % dtheta/dx
                dtetadx=calculotetadx(x,y);
                %==========================================================
                var1=((alfa*(r^alfa))/r^2)*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                var2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1));
                %
                part1=var1*x;
                part2=var2*dtetadx;
                % dtheta/dy
                dtetady=calculotetady(x,y);
                part3=var1*y;
                part4=var2*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=-0.7604;
                b3=-0.6495;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            %Obtain the flow rate
            F(ianalit,1) = dot(V,nij');
        end  %End of FOR (boundary edges)
        
        for ianalit=1:size(inedge,1)
            lef=inedge(ianalit,3);
            
            IJ=coord(inedge(ianalit,2),:)-coord(inedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(ianalit,2),:)+coord(inedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            
            teta1=calculoteta(x,y);
            alfa=0.6142;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-1.0546;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-0.4275;
                b2=0.2142;
                k11=10;
                k22=1000;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=-0.7604;
                b3=-0.6495;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            
            F(ianalit+size(bedge,1),1) = dot(V,nij');
        end
        
        vel=F;
        
        
    case 'zhangkobaise3'
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            teta=calculoteta(x,y);
            alfa=0.8866;
            r=sqrt(x^2+y^2);
            coef=alfa*(alfa-1)*(x^2+y^2)^(0.5*alfa-2);
            if elem(i,5)==1
                a1=1.0;
                b1=-0.3706;
                
                u(i,1)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
                kmap(1,:) = [1 1 0 0 1];
                fonte(i,1) =0;
            elseif elem(i,5)==2
                a2=-0.0144;
                b2=0.0022;
                
                u(i,1)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                kmap(2,:) = [2 10 0 0 100000];
                %=========================================================
                sum1= (x^2)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                sum2=(x*y)*(-2*b2*cos(alfa*teta)+2*a2*sin(alfa*teta));
                sum3=(y^2)*(-a2*cos(alfa*teta)-b2*sin(alfa*teta));
                sumtotal=sum1+sum2+sum3;
                fonte(i,1) =(10-100000)*coef*sumtotal*elemarea(i,1);
            else
                a3=0.7544;
                b3=-0.6564;
                
                u(i,1)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
                kmap(3,:) = [3 100 0 0 100];
                fonte(i,1) =0;
            end
        end  %End of FOR
        K=kmap;
        %% velocities
        for ianalit = 1:size(bedge,1)
            lef=bedge(ianalit,3);
            IJ=coord(bedge(ianalit,2),:)-coord(bedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            
            auxpoint=(coord(bedge(ianalit,2),:)+coord(bedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            teta1=calculoteta(x,y);
            alfa=0.8866;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-0.3706;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-0.0144;
                b2=0.0022;
                k11=10;
                k22=100000;
                % dtheta/dx
                dtetadx=calculotetadx(x,y);
                %==========================================================
                var1=((alfa*(r^alfa))/r^2)*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                var2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1));
                %
                part1=var1*x;
                part2=var2*dtetadx;
                % dtheta/dy
                dtetady=calculotetady(x,y);
                part3=var1*y;
                part4=var2*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=0.7544;
                b3=-0.6564;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            %Obtain the flow rate
            F(ianalit,1) = dot(V,nij');
        end  %End of FOR (boundary edges)
        
        for ianalit=1:size(inedge,1)
            lef=inedge(ianalit,3);
            
            IJ=coord(inedge(ianalit,2),:)-coord(inedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(ianalit,2),:)+coord(inedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            
            teta1=calculoteta(x,y);
            alfa=0.8866;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-0.3706;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-0.0144;
                b2=0.0022;
                k11=10;
                k22=100000;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=0.7544;
                b3=-0.6564;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            
            F(ianalit+size(bedge,1),1) = dot(V,nij');
        end
        
        vel=F;
        
    case 'zigzagfract'
        for i = 1:size(centelem,1)
            
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            if x<0.292
                if (0.15<x && x<0.292) && (0.48<y && y<0.52)
                    elem(i,5)=4;
                    
                    kmap(4,:) = [4 0.75025e-3 0.43258e-3 0.43258e-3 0.25025e-3];
                    %fonte(i,1)=100;
                else
                    elem(i,5)=1;
                    
                    kmap(1,:) = [1 0.75025e-3 0.43258e-3 0.43258e-3 0.25025e-3];
                    %fonte(i,1)=0;
                end
            elseif 0.292<x && x<0.711
                if (0.292<x && x<0.711) && (0.48<y && y<0.52)
                    elem(i,5)=5;
                    
                    kmap(5,:) = [5 0.25025e-3 -0.43258e-3 -0.43258e-3 0.75025e-3];
                    %fonte(i,1)=100;
                else
                    elem(i,5)=2;
                    
                    kmap(2,:) = [2 0.25025e-3 -0.43258e-3 -0.43258e-3 0.75025e-3];
                    %fonte(i,1)=0;
                end
            else
                
                if ((0.711<x && x<0.85) && (0.48<y && y<0.52))
                    elem(i,5)=6;
                    
                    kmap(6,:) = [6 0.25025e-3 0.43258e-3 0.43258e-3 0.75025e-3];
                    %fonte(i,1)=100;
                else
                    elem(i,5)=3;
                    
                    kmap(3,:) = [3 0.25025e-3 0.43258e-3 0.43258e-3 0.75025e-3];
                    %fonte(i,1)=0;
                end
            end
            
        end  %End of FOR
        K=kmap;
        
    case 'herbinhubert'
        
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            elem(i,5)=1;
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = 1.5;
            k(1,2) = 0.5;
            k(2,1) = 0.5;
            k(2,2) = 1.5;
            
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            u(i,1)= 16*x*(1-x)*y*(1-y);
            fonte(i,1) =1.5*32*(y*(1-y)+x*(1-x))-16*(1-2*x)*(1-2*y);
        end  %End of FOR
        K=kmap;
        
    case 'herbin'
        R=[0 -1 0; 1 0 0; 0 0 0];
        %Initialize a parameters
        alfa = 1000;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            elem(i,5)=i;
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            kaux(1,1) = (alfa*(x^2) + (y^2));
            kaux(1,2) = (alfa - 1)*(x*y);
            kaux(2,1) = (alfa - 1)*(x*y);
            kaux(2,2) = (alfa*(y^2) + (x^2));
            
            k = (1/((x^2) + (y^2))).*kaux;
            
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            u(i,1)= sin(pi*x)*sin(pi*y);
            fonte(i,1) =(pi/((centelem(i,1)^2) + ...
                (centelem(i,2)^2)))*(-(alfa - 1)*centelem(i,1)*...
                cos(pi*centelem(i,1))*(2*pi*centelem(i,2)*...
                cos(pi*centelem(i,2)) + sin(pi*centelem(i,2))) + ...
                sin(pi*centelem(i,1))*(-(alfa - 1)*centelem(i,2)*...
                cos(pi*centelem(i,2)) + (1 + alfa)*pi*...
                ((centelem(i,1)^2) + (centelem(i,2)^2))*...
                sin(pi*centelem(i,2))));
            
        end  %End of FOR
        K=kmap;
        %Initialize "velanalit"
        F = zeros(size(bedge,1) + size(inedge,1),1);
        %Calculate the analitical velocity
        %Boundary edges
        for ianalit = 1:size(bedge,1)
            
            %Obtain the vector velocity (-KnablaP)
            %     V = [(-pi/((overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2)))*...
            %         ((delta - 1)*overedgecoord(ianal,1)*overedgecoord(ianal,2)*...
            %         cos(pi*overedgecoord(ianal,2))*sin(pi*overedgecoord(ianal,1)) + ...
            %         (delta*(overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2))*...
            %         cos(pi*overedgecoord(ianal,1))*sin(pi*overedgecoord(ianal,2)));
            %         (-pi/((overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2)))*...
            %         ((delta - 1)*overedgecoord(ianal,1)*overedgecoord(ianal,2)*...
            %         cos(pi*overedgecoord(ianal,1))*sin(pi*overedgecoord(ianal,2)) + ...
            %         (delta*(overedgecoord(ianal,2)^2) + (overedgecoord(ianal,1)^2))*...
            %         cos(pi*overedgecoord(ianal,2))*sin(pi*overedgecoord(ianal,1)));
            %         0];
            IJ=coord(bedge(ianalit,2),:)-coord(bedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            
            auxpoint=(coord(bedge(ianalit,2),:)+coord(bedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            dpdx=pi*cos(pi*x)*sin(pi*y);
            dpdy=pi*sin(pi*x)*cos(pi*y);
            k11=(alfa*x^2 + y^2)/(x^2+y^2);
            k12= ((alfa-1)*x*y)/(x^2+y^2);
            k22=(x^2+alfa*y^2)/(x^2+y^2);
            
            V=-[k11*dpdx+k12*dpdy k12*dpdx+k22*dpdy 0];
            %Obtain the flow rate
            F(ianalit) = dot(V,nij');
        end  %End of FOR (boundary edges)
        
        %Internal edges
        for ianalit = 1:size(inedge,1)
            IJ=coord(inedge(ianalit,2),:)-coord(inedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(ianalit,2),:)+coord(inedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            dpdx=pi*cos(pi*x)*sin(pi*y);
            dpdy=pi*sin(pi*x)*cos(pi*y);
            k11=(alfa*x^2 + y^2)/(x^2+y^2);
            k12= ((alfa-1)*x*y)/(x^2+y^2);
            k22=(x^2+alfa*y^2)/(x^2+y^2);
            
            V=-[k11*dpdx+k12*dpdy k12*dpdx+k22*dpdy 0];
            %Obtain the flow rate
            F(size(bedge,1) + ianalit) = ...
                dot(V,nij');
        end
        vel=F;
    case 'gaowu6'
        fonte=0;
        elem(:,5)=1;
        u=0;
        %Initialize "R":
        R = zeros(2);
        % problema Gao e Wu 2013
        k = [1 0; 0 1e-3];
        % problema TEREKHOV
        %k = [50 0; 0 1];
        %Fill "R"
        R(1,1) = cosd(67.5);
        R(1,2) = sind(67.5);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        K=kmap;
    case 'benchmar5_6'
        %% problem 5.6 CASO 2
        u=0;
        for ielem=1:size(elem,1)
            if (0<centelem(ielem,1) && centelem(ielem,1)<0.5) && ...
                    (0<centelem(ielem,2) && centelem(ielem,2)<0.5)
                elem(ielem,5)=ielem;
                theta=pi/6;
                k=[1000 0;0 1];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0.5<centelem(ielem,1) && centelem(ielem,1)<1) && ...
                    (0<centelem(ielem,2) && centelem(ielem,2)<0.5)
                elem(ielem,5)=ielem;
                theta=-pi/6;
                k=[10 0;0 1];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0.5<centelem(ielem,1) && centelem(ielem,1)<1) && ...
                    (0.5<centelem(ielem,2) && centelem(ielem,2)<1)
                elem(ielem,5)=ielem;
                theta=pi/6;
                k=[1000 0;0 1];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0<centelem(ielem,1) && centelem(ielem,1)<0.5) && ...
                    (0.5<centelem(ielem,2) && centelem(ielem,2)<1)
                elem(ielem,5)=ielem;
                theta=-pi/6;
                k=[10 0;0 1];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            end
            if (7/18<centelem(ielem,1) || 7/18==centelem(ielem,1))&& ...
                    (centelem(ielem,1)<11/18 || 11/18==centelem(ielem,1))...
                    && (7/18<centelem(ielem,2) || centelem(ielem,2)==7/18) && ...
                    (centelem(ielem,2)<11/18 || centelem(ielem,2)==11/18)
                fonte(ielem,1)=81/4;
            else
                fonte(ielem,1)=0;
            end
            K=kmap;
            
        end
        fonte=fonte.*elemarea;
    case 'benchmar5_7' %
        %% problem 5.6 CASO 1
        u=0;
        k=[1000 0;0 1];
        for ielem=1:size(elem,1)
            if (0<centelem(ielem,1) && centelem(ielem,1)<0.5) && ...
                    (0<centelem(ielem,2) && centelem(ielem,2)<0.5)
                elem(ielem,5)=ielem;
                theta=pi/6;
                
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0.5<centelem(ielem,1) && centelem(ielem,1)<1) && ...
                    (0<centelem(ielem,2) && centelem(ielem,2)<0.5)
                elem(ielem,5)=ielem;
                theta=-pi/6;
                
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0.5<centelem(ielem,1) && centelem(ielem,1)<1) && ...
                    (0.5<centelem(ielem,2) && centelem(ielem,2)<1)
                elem(ielem,5)=ielem;
                theta=pi/6;
                
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0<centelem(ielem,1) && centelem(ielem,1)<0.5) && ...
                    (0.5<centelem(ielem,2) && centelem(ielem,2)<1)
                elem(ielem,5)=ielem;
                theta=-pi/6;
                
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            end
            if (7/18<centelem(ielem,1) || 7/18==centelem(ielem,1))&& ...
                    (centelem(ielem,1)<11/18 || 11/18==centelem(ielem,1))...
                    && (7/18<centelem(ielem,2) || centelem(ielem,2)==7/18) && ...
                    (centelem(ielem,2)<11/18 || centelem(ielem,2)==11/18)
                fonte(ielem,1)=81/4;
            else
                fonte(ielem,1)=0;
            end
        end
        K=kmap;
        fonte=fonte.*elemarea;
    case 'edwards'
        %% This problem is adapted of (Rogers and Edwards 1998)
        a1=1;
        f=((4*a1)/(((1/50)-2)*(1/10)+1));
        b2=((1/10)-1)*f;
        c2=f;
        d2=-c2*(1/10);
        c1=(1/50)*(1/10)*f;
        d1=d2;
        % calculo da solução exata, termo fonte e adequa o tensor de
        % permeabilidade
        for ielem=1:size(elem,1)
            x=centelem(ielem,1);
            y=centelem(ielem,2);
            if centelem(ielem,1)<0.5
                u(ielem,1)=c1*x^2+d1*y^2;
                kmap(ielem,1:5)=[ielem 50 0 0 1];
                elem(ielem,5)=ielem;
                fonte(ielem,1)=-100*c1-2*d1;
                normKmap(ielem,1)=norm([50 0; 0 1]);
            else
                u(ielem,1)=1+b2*x+c2*x^2+d2*y^2;
                kmap(ielem,1:5)=[ielem 1 0 0 10];
                elem(ielem,5)=ielem;
                fonte(ielem,1)=-2*c2-20*d2;
                normKmap(ielem,1)=norm([1 0; 0 10]);
            end
        end
        % calculo das velocidades
        R=[0 -1 0; 1 0 0; 0 0 0];
        for i=1:size(bedge,1)
            
            IJ=coord(bedge(i,2),:)-coord(bedge(i,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(bedge(i,2),:)+coord(bedge(i,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            if x<0.5
                a= [-100*c1*x -2*d1*y 0];
            else
                a= [-b2-2*c2*x -20*d2*y 0];
            end
            F(i,1)=dot(a,nij');
            
        end
        
        for i=1:size(inedge,1)
            
            IJ=coord(inedge(i,2),:)-coord(inedge(i,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(i,2),:)+coord(inedge(i,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            if x<0.5
                a= [-100*c1*x -2*d1*y 0];
            else
                a= [-b2-2*c2*x -20*d2*y 0];
            end
            F(i+size(bedge,1),1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'FriisEdwards'
        u=0;
        fonte=zeros(size(elem,1),1);
        %[element1]=elementsneigbormiddleface;
        % 253   five_spotMALHAV1
        % 5178  five_spotMALHAV2
        % 21517 five_spotMALHAV3
        % 20064 five_spotMALHAV4
        % 9985 MESH_06TRIANGULARDISTORTED_96PROBLEM431
        
        fonte(9985,1)=248000000;
        
        k = [3000 0; 0 1];
        %Fill "R"
        theta=-pi/6;
        R(1,1) = cos(theta);
        R(1,2) = sin(theta);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        elem(:,5)=1;
        K=kmap;
    case 'shenyuan16'
        
        for i=1:size(elem,1)
            x= centelem(i,1);
            y=centelem(i,2);
            k = [1+2*x^2+y^2 0; 0 1+x^2+2*y^2];
            %Fill "R"
            theta=5*pi/12;
            R(1,1) = cos(theta);
            R(1,2) = sin(theta);
            R(2,1) = -R(1,2);
            R(2,2) = R(1,1);
            %Fill "k" turning the tensor
            A=inv(R);
            k = A*k*R;
            %Buld "kmap" again
            kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            % valido
            t1=2.46711*(x^2-y^2)*cos(pi*x)*cos(pi*y)-pi^2*(1+1.0669*x^2+1.93289*y^2)*sin(pi*x)*sin(pi*y)+...
                1.57061*x*sin(pi*x)*cos(pi*y)+6.70353*x*cos(pi*x)*sin(pi*y);
            
            t2=2.46711*(x^2-y^2)*cos(pi*x)*cos(pi*y)-pi^2*(1+1.93289*x^2+1.0669*y^2)*sin(pi*x)*sin(pi*y)+...
                6.70353*y*sin(pi*x)*cos(pi*y)-1.57061*y*cos(pi*x)*sin(pi*y);
            fonte(i,1)=-(t1+t2)*elemarea(i,1) ;
            
            u(i,1)=sin(pi*x)*sin(pi*y);
            elem(i,5)=i;
        end
        for iface=1:size(bedge,1)+size(inedge,1)
            R1=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                t0= pi*(1+1.0670*x^2+1.9330*y^2)*cos(pi*x)*sin(pi*y)+pi*(x^2-y^2)*0.25*sin(pi*x)*cos(pi*y);
                t01=pi*(1+1.9330*x^2+1.0670*y^2)*sin(pi*x)*cos(pi*y)+pi*(x^2-y^2)*0.25*cos(pi*x)*sin(pi*y);
                
                % anterior
                
                %t0= pi*(1+1.0669*x^2+1.93289*y^2)*cos(pi*x)*sin(pi*y)+pi*(x^2-y^2)*0.24997*sin(pi*x)*cos(pi*y);
                %t01=pi*(1+1.93289*x^2+1.0669*y^2)*sin(pi*x)*cos(pi*y)+pi*(x^2-y^2)*0.24997*cos(pi*x)*sin(pi*y);
                
                
                a=[-t0, -t01, 0];
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                t0= pi*(1+1.0670*x^2+1.9330*y^2)*cos(pi*x)*sin(pi*y)+pi*(x^2-y^2)*0.25*sin(pi*x)*cos(pi*y);
                t01=pi*(1+1.9330*x^2+1.0670*y^2)*sin(pi*x)*cos(pi*y)+pi*(x^2-y^2)*0.25*cos(pi*x)*sin(pi*y);
                
                % anterior
                %t0= pi*(1+1.0669*x^2+1.93289*y^2)*cos(pi*x)*sin(pi*y)+pi*(x^2-y^2)*0.24997*sin(pi*x)*cos(pi*y);
                %t01=pi*(1+1.93289*x^2+1.0669*y^2)*sin(pi*x)*cos(pi*y)+pi*(x^2-y^2)*0.24997*cos(pi*x)*sin(pi*y);
                
                
                a=[-t0, -t01, 0];
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'gaowu5'
        fonte=0;
        elem(:,5)=1;
        u=0;
        %Initialize "R":
        R = zeros(2);
        k = [1 0; 0 1e-4];
        %Fill "R"
        R(1,1) = cosd(40);
        R(1,2) = sind(40);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        K=kmap;
        
    case 'edqueiroz'
        % problema artigo quiroz et al.
        epsilon=1e3;
        for ielem=1:size(elem,1)
            if centelem(ielem,1)<0.5 ||  centelem(ielem,1)==0.5
                theta=0.5*pi;
                k=[100 0;0 0.01];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
                elem(ielem,5)=ielem;
            else
                x1= centelem(ielem,1)+1e-3;
                y1= centelem(ielem,2)+1e-3;
                
                kmap(ielem,1:5)=[ielem y1^2+epsilon*x1^2 -(1-epsilon)*x1*y1 -(1-epsilon)*x1*y1 x1^2+epsilon*y1^2 ];
                elem(ielem,5)=ielem;
                
            end
        end
        K=kmap;
        fonte=0;
        u=0;
        %% use se deseja resolver o problema de Queiroz et al 2014
        % unstructured mesh com furo reordenando o sentido da fronterira no
        % contorno interior
        % unstructured mesh com furo reordenando o sentido da fronterira no
        % contorno interior
%          x=bedge(71:78,1); % UTILIZE Benchmark23_3_18_18.msh
%          y=bedge(71:78,2);
%          bedge(71:78,1)=y;
%          bedge(71:78,2)=x;
%          bedge(71:78,4:5)=102; % 18x18
%          bcflag(2,1)=102;
%          bcflag(2,2)=2;
        %=====================================
%          x=bedge(73:80,1); % UTILIZE MeshTri18sym.msh
%          y=bedge(73:80,2);
%          bedge(73:80,1)=y;
%          bedge(73:80,2)=x;
%          bedge(73:80,4:5)=102; % 18x18
%          bcflag(2,1)=102;
%          bcflag(2,2)=2;
        
        %=====================================
        %
        x=bedge(135:150,1); % UTILIZE Benchmark23_3_FINA36.msh
        y=bedge(135:150,2);
        bedge(135:150,1)=y;
        bedge(135:150,2)=x;
        bedge(135:150,4:5)=102; % 36x36
        bcflag(2,1)=102;
        bcflag(2,2)=2;
        %===============================================================
        %  x=bedge(145:160,1); % UTILIZE MeshTri36sym.msh
        %  y=bedge(145:160,2);
        %  bedge(145:160,1)=y;
        %  bedge(145:160,2)=x;
        %  bedge(145:160,4:5)=102; % 36x36
        %  bcflag(2,1)=102;
        %  bcflag(2,2)=2;
        
        %===============================================================
        % unstructured mesh com furo reordenando o sentido da fronterira no
        % contorno interior
%         x=bedge(289:320,1);
%         y=bedge(289:320,2);
%         bedge(289:320,1)=y;
%         bedge(289:320,2)=x;
%         bedge(289:320,4:5)=102; % benchmark23_3 72x72
%         bcflag(2,1)=102;
%         bcflag(2,2)=2;
        
        % unstructured mesh com furo reordenando o sentido da fronterira no
        % contorno interior
        %  x=bedge(577:640,1);
        %  y=bedge(577:640,2);
        %  bedge(577:640,1)=y;
        %  bedge(577:640,2)=x;
        %  bedge(577:640,4:5)=102; % benchmark23_3 144x144
        %  bcflag(2,1)=102;
        %  bcflag(2,2)=2;
        
    case 'homogeneo'
        %% meio isotropico homogeneo
        for i=1:size(centelem,1)
            
            kmap(1,1:5)=[1 1 0 0 1];
            elem(i,5)=1;
            
            u(i,1)=1-centelem(i,1);
            fonte=0;
        end
        
        for iface=1:size(bedge,1)+size(inedge,1)
            R1=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                
                a=[1, 0, 0];
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                
                a=[1, 0, 0];
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'heterogeneo'
        for i=1:size(elem,1)
            if centelem(i,1)<0.5
                kmap(1,1:5)=[1 1 0 0 1];
                elem(i,5)=1;
                u(i,1)=(4/3)*centelem(i,1);
            else
                
                kmap(2,1:5)=[1 2 0 0 2];
                
                elem(i,5)=2;
                u(i,1)=(2/3)*centelem(i,1)+(1/3);
            end
            fonte=0;
        end
        for iface=1:size(bedge,1)+size(inedge,1)
            R1=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                if IJ(1,1)<0.5
                    
                    a=[-4/3, 0, 0];
                else
                    a=[-4/3, 0, 0];
                end
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                
                if IJ(1,1)<0.5
                    
                    a=[-4/3, 0, 0];
                else
                    a=[-4/3, 0, 0];
                end
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'crumpton'
        
        alfa=1000;
        u=zeros(size(elem,1),1);
        for ielem=1:size(elem,1)
            if centelem(ielem,1)<0 || centelem(ielem,1)==0
                u(ielem)=(2*sin(centelem(ielem,2))+cos(centelem(ielem,2)))*alfa*centelem(ielem,1)+sin(centelem(ielem,2))+3*alfa;
                kmap(1,1:5)=[1 1 0 0 1];
                elem(ielem,5)=1;
                fonte(ielem,1)=(2*sin(centelem(ielem,2))+cos(centelem(ielem,2)))*alfa*centelem(ielem,1)+sin(centelem(ielem,2));
            else
                u(ielem)=exp(centelem(ielem,1))*sin(centelem(ielem,2))+3*alfa;
                kmap(2,1:5)=[2 alfa*2 alfa*1 alfa*1 alfa*2];
                elem(ielem,5)=2;
                fonte(ielem,1)=-2*alfa*exp(centelem(ielem,1))*cos(centelem(ielem,2));
            end
        end
        
        % calculo dos fluxos analiticos
        R=[0 -1 0; 1 0 0; 0 0 0];
        for i=1:size(bedge,1)
            
            IJ=coord(bedge(i,2),:)-coord(bedge(i,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(bedge(i,2),:)+coord(bedge(i,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            if x<0 || x==0
                a= [-(2*sin(y)+ cos(y))*alfa -(2*cos(y)- sin(y))*alfa*x-cos(y) 0];
            else
                a= alfa*exp(x)*[-(2*sin(y)+ cos(y)) -(sin(y)+ 2*cos(y)) 0];
            end
            F(i,1)=dot(a,nij');
            
        end
        
        for i=1:size(inedge,1)
            
            IJ=coord(inedge(i,2),:)-coord(inedge(i,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(i,2),:)+coord(inedge(i,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            if x<0 || x==0
                a= [-(2*sin(y)+ cos(y))*alfa -(2*cos(y)- sin(y))*alfa*x-cos(y) 0];
            else
                a= alfa*exp(x)*[-(2*sin(y)+ cos(y)) -(sin(y)+ 2*cos(y)) 0];
            end
            F(i+size(bedge,1),1)=dot(a,nij');
            
        end
        K=kmap;
        vel=F;
    case 'gaowu1'
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            if x<0.5 || x==0.5
                fonte(i,1)=-exp(-62.8319*((-0.5+y)^2))*(127304-782890*y+(1.55068*(10^6))*...
                    (y^2)-992201*(y^3)+x*(-222452+(1.43303*10^6)*y-(2.96871*10^6)*(y^2)+...
                    (1.9844*10^6)*(y^3)));
                % obtido no wolfram
                
                %                fonte(i,1)=-exp(-62.8319*((y-0.5)^2))*(128821-789183*y+(1.557*(10^6))*...
                %                     (y^2)-992204*(y^3)+x*(-222452+(1.43303*10^6)*y-(2.96871*10^6)*(y^2)+...
                %                     (1.9844*10^6)*(y^3)));
                
                u(i,1)=(1+(x-0.5)*(0.1+8*pi*(y-0.5)))*exp(-20*pi*((y-0.5)^2));
                elem(i,5)=1;
                kmap(1,:)=[1 10 2 2 5];
            else
                fonte(i,1)=(-exp(x-62.8319*((y-0.5)^2))*(2318.87-9577.95*y+9577.95*(y^2)));
                % obtido por wolfram
                %fonte(i,1)=-exp(x-62.8319*((y-0.5)^2))*(23.789231-95.7796*y+95.7796*(y^2));
                u(i,1)=exp(x-0.5)*exp(-20*pi*((y-0.5)^2));
                
                elem(i,5)=2;
                kmap(2,:)=[2 1 0 0 1];
            end
        end
        K=kmap;
    case 'gaowu2'
        c=centelem;
        R=[0 -1 0; 1 0 0; 0 0 0];
        for i=1:size(elem,1)
            x=c(i,1);
            y=c(i,2);
            % valido
            fonte(i,1)= (-0.25*(-3*(-1 + c(i,1))^2*(-1 + c(i,2)) + cos((-1 + c(i,1))*(-1 + c(i,2)))* csc(1)) -  0.25*(-1 + c(i,2))*(-3*(-1 + c(i,1))^2 - (-1 + c(i,1))*csc(1)*sin((-1 + c(i,1))*(-1 + c(i,2)))) - ...
                0.75*(-2*(-1 + c(i,1))^3 - (-1 + c(i,1))^2*csc(1)*sin((-1 + c(i,1))*(-1 + c(i,2)))) - 0.75*(-1 + c(i,2))*(-6*(-1 + c(i,1))*(-1 + c(i,2)) - (-1 + c(i,2))*csc(1)* sin((-1 + c(i,1))*(-1 + c(i,2)))) - ...
                0.25*(-6*(-1 + c(i,1))^2 *(-1 + c(i,2)) +cos((-1 + c(i,1))*(-1 + c(i,2)))*csc(1) - (-1 + c(i,1))*(-1 + c(i,2))*csc(1)*sin((-1 + c(i,1))*(-1 + c(i,2)))));
            %             fonte(i,1)=-0.25*(-3*((x - 1)^2)*(y - 1) + ...
            %                     cos((x - 1)*(y - 1))*csc(1)) - ...
            %                     0.25*(y - 1)*(-3*((x - 1)^2) - ...
            %                     (x - 1)*csc(1)*sin((x - 1)*(y - 1))) - 0.75*(-2*((x - 1)^3) - ...
            %                     ((x - 1)^2)*csc(1)*sin((x - 1)*(y - 1))) - ...
            %                     0.75*(y - 1)*(-6*(x - 1)*(y - 1) - ...
            %                     (y - 1)*csc(1)*sin((x - 1)*(y - 1))) - ...
            %                     0.25*(-6*((x - 1)^2)*(y - 1) + cos((x - 1)*(y - 1))*csc(1) - ...
            %                     (x - 1)*(y - 1)*csc(1)*sin((x - 1)*(y - 1)));
            
            x= centelem(i,1);
            y=centelem(i,2);
            u(i,1)=0.5*((sin((1-x)*(1-y))/(sin(1)))+(1-x)^3*(1-y)^2);
            elem(i,5)=1;
        end
        
        for iface=1:size(bedge,1)+size(inedge,1)
            
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
                dpdx=((1-y)/2)*(((-cos((1-x)*(1-y)))/sin(1))-3*(1-x)^2*(1-y));
                dpdy=((1-x)/2)*(((-cos((1-x)*(1-y)))/sin(1))-2*(1-x)^2*(1-y));
                a=-[1.5*dpdx+0.5*dpdy 0.5*dpdx+1.5*dpdy 0];
                
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
                
                dpdx=((1-y)/2)*(((-cos((1-x)*(1-y)))/sin(1))-3*(1-x)^2*(1-y));
                dpdy=((1-x)/2)*(((-cos((1-x)*(1-y)))/sin(1))-2*(1-x)^2*(1-y));
                a=-[1.5*dpdx+0.5*dpdy 0.5*dpdx+1.5*dpdy 0];
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        
        kmap=[1 1.5 0.5 0.5 1.5];
        K=kmap;
        vel=F;
    case 'crumptonhyman'
        
        for i=1:size(elem,1)
            x= centelem(i,1);
            y=centelem(i,2);
            % valido
            fonte(i,1)= -2*(1+x^2+x*y+y^2)*exp(x*y);
            u(i,1)=exp(x*y);
            elem(i,5)=1;
        end
        kmap=[1 2 1 1 2];
    case 'gaowu3'
        alfa=1000;
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            r=exp(-20*pi*((x-0.5)^2+(y-0.5)^2));
            %% a mano
            %             fonte(i,1)=-r*(-40*pi*(3*alfa-1)*(x*(x-0.5)+y*(y-0.5))+...
            %                 3200*pi^2*x*y*(x-0.5)*(y-0.5)*(alfa-1)+...
            %                 (alfa*x^2+y^2)*(-40*pi+1600*pi^2*(x-0.5)^2)+...
            %                 (x^2+alfa*y^2)*(-40*pi+1600*pi^2*(y-0.5)^2));
            %% wolfram alpha
            
            parte1=(y^2+alfa*x^2)*r*(- 125.664 + 15791.4*(x-0.5)^2)+...
                r*15791.4*(alfa-1)*(y-0.5)*y*x*(x-0.5)- 251.328*alfa*x*(x-0.5)*r-...
                125.664*(alfa-1)*(y-0.5)*y*r ;
            
            parte2=r*(x^2+alfa*y^2)*(-125.664 + 15791.4*(y-0.5)^2)+...
                r*15791.4*(alfa-1)*(x-0.5)*x*y*(y-0.5) - 251.328*alfa*y*(y-0.5)*r-...
                125.664*(alfa-1)*(x-0.5)*x*r;
            
            fonte(i,1)=-parte1-parte2;
            u(i,1)=exp(-20*pi*((x-0.5)^2 + (y-0.5)^2))+1e-3;
            elem(i,5)=i;
            kmap(i,1:5)=[i alfa*x^2+y^2 (alfa-1)*x*y (alfa-1)*x*y x^2+alfa*y^2];
        end
        
    case 'lepotier'
        R=[0 1 0; -1 0 0; 0 0 0];
        alfa=1e-3;
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            uu=sin(pi*x)*sin(pi*y);
            r=cos(pi*x)*sin(pi*y);
            s=sin(pi*x)*cos(pi*y);
            %% no artigo de Le potier
            %fonte(i,1)=elemarea(i,1)*uu*((1+alfa)*pi^2*(x^2+y^2))-r*((1-3*alfa)*pi*x)-s*((1-3*alfa)*pi*y)-...
            %    cos(pi*x)*cos(pi*y)*(2*pi^2*(1-alfa)*x*y);
            
            %k=1;
            fonte(i,1)=elemarea(i,1)*(pi/(x^2+y^2))*(-(alfa - 1)*x*...
                cos(pi*x)*(2*pi*y*...
                cos(pi*y) + sin(pi*y)) + ...
                sin(pi*x)*(-(alfa - 1)*y*...
                cos(pi*y) + (1 + alfa)*pi*((x^2) + (y^2))*...
                sin(pi*y)));
            
            k = (1/((x^2) + (y^2)));
            
            % para artigo de Le-potier considere
            %u(i,1)=sin(pi*x)*sin(pi*y);
            % para artigo de zhang Kobaise 2019 considere:
            u(i,1)=sin(pi*x)*sin(pi*y)+1;
            elem(i,5)=i;
            kmap(i,1:5)=[i k*(alfa*x^2+y^2) k*((alfa-1)*x*y) k*((alfa-1)*x*y) k*(x^2+alfa*y^2)];
        end
        K=kmap;
        for iface=1:size(bedge,1)+size(inedge,1)
            
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
            end
            
            dpdx=pi*cos(pi*x)*sin(pi*y);
            dpdy=pi*sin(pi*x)*cos(pi*y);
            k11=(y^2+alfa*x^2)/(x^2 + y^2);
            k12=((alfa-1)*x*y)/(x^2 + y^2);
            k22=(alfa*y^2+x^2)/(x^2 + y^2);
            
            a=-[k11*dpdx+k12*dpdy  k12*dpdx+k22*dpdy 0];
            
            F(iface,1)=dot(a,(R*IJ')');
            
        end
        vel=F;
    case 'gaowu4'
        alfa=1e-6;
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            uu=sin(pi*x)*sin(pi*y);
            r=cos(pi*x)*sin(pi*y);
            s=sin(pi*x)*cos(pi*y);
            d=x^2+y^2;
            %% obtido por wolfram alpha1
            %             part1= (1/d)*(-pi^2*uu*(alfa*x^2+y^2)+ pi^2*(alfa-1)*x*y*cos(pi*x)*cos(pi*y)+ ...
            %                 pi*(alfa-1)*y*s+ 2*pi*alfa*x*r)- (1/d^2)*(2*x*(pi*r*(alfa*x^2+y^2)+pi*(alfa-1)*x*y*s));
            %
            %             part2=(1/d)*(-pi^2*uu*(alfa*y^2+x^2)+ pi^2*(alfa-1)*x*y*cos(pi*x)*cos(pi*y)+ ...
            %                 2*pi*alfa*y*s+ pi*(alfa-1)*x*r)-(1/d^2)*(2*y*(pi*s*(alfa*y^2+x^2)+pi*(alfa-1)*x*y*r));
            %% wolfram alpha2
            
            part1= (1/d)*(-pi^2*uu*(alfa*x^2+y^2)- pi^2*(alfa-1)*x*y*uu+ ...
                pi*(alfa-1)*x*s+ 2*pi*alfa*x*r)- (1/d^2)*(2*x*pi*(r*(alfa*x^2+y^2)+(alfa-1)*y^2*s));
            
            
            
            part2=(1/d)*(-pi^2*uu*(alfa*y^2+x^2)- pi^2*(alfa-1)*x*y*uu+ ...
                2*pi*alfa*y*s+ pi*(alfa-1)*y*r) -   (1/d^2)*(2*y*pi*(s*(alfa*y^2+x^2)+(alfa-1)*x^2*r));
            
            fonte(i,1)=-part1-part2;
            
            u(i,1)=sin(pi*x)*sin(pi*y);
            elem(i,5)=i;
            kmap(i,1:5)=[i (alfa*x^2+y^2)/d ((alfa-1)*x*y)/d ((alfa-1)*x*y)/d (x^2+alfa*y^2)/d];
            
        end
    case 'lipnikov1'
        
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            if x<0.5 || x==0.5
                fonte(i,1)=4;
                u(i,1)=1-2*y^2+4*x*y+6*x+2*y;
                elem(i,5)=1;
                kmap(1,1:5)=[1 1 0 0 1];
            else
                fonte(i,1)=-5.6 ;
                u(i,1)=-2*y^2+1.6*x*y-0.6*x+3.2*y+4.3;
                elem(i,5)=2;
                kmap(2,1:5)=[2 10 3 3 1];
            end
            
        end
        
    case 'lipnikov2'
        u=0;
        fonte=0;
        k = [1000 0; 0 1];
        %Fill "R"
        %theta=5*pi/6;
        theta=67.5;
        R(1,1) = cosd(theta);
        R(1,2) = sind(theta);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        elem(:,5)=1;
        K=kmap;
        
    case 'lipnikov3'
        u=0;
        fonte=0;
        k = [100000 0; 0 1];
        %Fill "R"
        %theta=pi/4;
        theta=-pi/6;
        R(1,1) = cos(theta);
        R(1,2) = sin(theta);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        elem(:,5)=1;
        K=kmap;
    case 'herbinhubert'
        %Initialize a parameter
        epsilon = 1e3;
        theta = 0.5*pi;
        k = [100 0; 0 0.01];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
            [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Define "x1" and "y1"
                x1 = x + 1e-3;
                y1 = y + 1e-3;
                %Definition of permeability components
                k(1,1) = ((y1^2) + epsilon*(x1^2));
                k(1,2) = -(1 - epsilon)*(x1*y1);
                k(2,1) = -(1 - epsilon)*(x1*y1);
                k(2,2) = ((x1^2) + epsilon*(y1^2));
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR
    case 'guangwei'
        u=0;
        alfa=5*1e-3;
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            kmap(i,1:5)=[i (y^2+alfa*x^2)+alfa -(1-alfa)*x*y -(1-alfa)*x*y  y^2*alfa+x^2+alfa];
            elem(i,5)=i;
            if (x>3/8 || x==3/8)&&(x<5/8 || x==5/8)&&(y>3/8 || y==3/8)&&(y<5/8 || y==5/8)
                fonte(i,1)=1;
            else
                fonte(i,1)=0;
            end
        end
        K=kmap;
    case 'guangwei1'
        
        u=0;
        %Initialize "R":
        R = zeros(2);
        k = [10000 0; 0 1];
        %Fill "R"
        R(1,1) = cos(pi/6);
        R(1,2) = sin(pi/6);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            elem(i,5)=1;
            %             if (7/18<centelem(i,1) || 7/18==centelem(i,1))&& ...
            %                     (centelem(i,1)<11/18 || 11/18==centelem(i,1))...
            %                     && (7/18<centelem(i,2) || centelem(i,2)==7/18) && ...
            %                     (centelem(i,2)<11/18 || centelem(i,2)==11/18)
            if (50/101<centelem(i,1) || 50/101==centelem(i,1))&& ...
                    (centelem(i,1)<51/101 || 51/101==centelem(i,1))...
                    && (50/101<centelem(i,2) || centelem(i,2)==50/101) && ...
                    (centelem(i,2)<51/101 || centelem(i,2)==51/101)
                fonte(i,1)=1021;
            else
                fonte(i,1)=0;
            end
        end
        
    case 'gaowu7'
        auxbedge1=bedge(:,1);
        auxbedge2=bedge(:,2);
        
        bedge(:,1)=auxbedge2;
        bedge(:,2)=auxbedge1;
        
        u=0;
        alfa=0.2;
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            phi1=y-alfa*(x-0.5)-0.475;
            phi2=phi1-0.05;
            theta=atand(alfa);
            if phi1<0 || phi1==0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [1 0; 0 0.1];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)=0;
            elseif phi1>0 && phi2<0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [100 0; 0 10];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)=0;
            elseif phi2>0 || phi2==0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [1 0; 0 0.1];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)= 0;
            end
            
            u(i,1)=-x-alfa*y;
            KK=kmap(i,2:5);
            r=[KK(1,1) KK(1,2);
                KK(1,3) KK(1,4)];
            normKmap(i,1)=norm(r);
        end
        K=kmap;
        
    case 'gaowu8'
        auxbedge1=bedge(:,1);
        auxbedge2=bedge(:,2);
        
        bedge(:,1)=auxbedge2;
        bedge(:,2)=auxbedge1;
        
        u=0;
        alfa=0.2;
        fonte=0;
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            phi1=y-alfa*(x-0.5)-0.475;
            phi2=phi1-0.05;
            
            if phi1<0
                elem(i,5)=i;
                %Buld "kmap" again
                kmap(i,1:5) = [i 1 0 0 1];
                u(i,1)=-phi1;
            elseif phi1>0 && phi2<0
                elem(i,5)=i;
                %Buld "kmap" again
                kmap(i,1:5) = [i 0.01 0 0 0.01];
                u(i,1)=-100*phi1;
                
            elseif phi2>0
                elem(i,5)=i;
                
                %Buld "kmap" again
                kmap(i,1:5) = [i 1 0 0 1];
                u(i,1)=-phi2-5;
            end
            KK=kmap(i,2:5);
            r=[KK(1,1) KK(1,2);
                KK(1,3) KK(1,4)];
            normKmap(i,1)=norm(r);
        end
        for iface=1:size(bedge,1)+size(inedge,1)
            R=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                phi1=y-alfa*(x-0.5)-0.475;
                phi2=phi1-0.05;
                if phi1<0
                    
                    %Initialize "R":
                    
                    k = [1 0; 0 1];
                    
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0;k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[alfa;-1;0];
                elseif phi1>0 && phi2<0
                    
                    %Initialize "R":
                    
                    k = [0.01 0; 0 0.01];
                    
                    %Buld "kmap" again
                    KK= [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[100*alfa;-100;0];
                elseif phi2>0
                    
                    %Initialize "R":
                    
                    k = [1 0; 0 1];
                    
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[alfa;-1;0];
                end
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
                phi1=y-alfa*(x-0.5)-0.475;
                phi2=phi1-0.05;
                theta=atand(alfa);
                if phi1<0
                    
                    %Initialize "R":
                    
                    k = [1 0; 0 1];
                    
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0;k(2,1) k(2,2) 0; 0 0 0];
                    
                elseif phi1>0 && phi2<0
                    
                    %Initialize "R":
                    
                    k = [0.01 0; 0 0.01];
                    
                    %Buld "kmap" again
                    KK= [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[100*alfa;-100;0];
                    
                elseif phi2>0
                    
                    %Initialize "R":
                    
                    k = [1 0; 0 1];
                    
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[alfa;-1;0];
                end
                
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'gaowu9'
        % adaptado de Gao e Wu, 2014
        auxbedge1=bedge(:,1);
        auxbedge2=bedge(:,2);
        
        bedge(:,1)=auxbedge2;
        bedge(:,2)=auxbedge1;
        
        u=0;
        alfa=0.2;
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            phi1=y-alfa*(x-0.5)-0.475;
            phi2=phi1-0.05;
            theta=atand(alfa);
            if phi1<0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [1 0; 0 0.1];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                
            elseif phi1>0 && phi2<0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [100 0; 0 10];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)=0;
            elseif phi2>0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [1 0; 0 0.1];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)= 0;
            end
            
            u(i,1)=2-x-alfa*y;
            KK=kmap(i,2:5);
            r=[KK(1,1) KK(1,2);
                KK(1,3) KK(1,4)];
            normKmap(i,1)=norm(r);
        end
        
        for iface=1:size(bedge,1)+size(inedge,1)
            R=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                phi1=y-alfa*(x-0.5)-0.475;
                phi2=phi1-0.05;
                theta=atand(alfa);
                if phi1<0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [1 0; 0 0.1];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0;k(2,1) k(2,2) 0; 0 0 0];
                    
                elseif phi1>0 && phi2<0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [100 0; 0 10];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK= [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    
                elseif phi2>0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [1 0; 0 0.1];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    
                end
                
                a=-KK*[-1;-0.2;0];
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
                phi1=y-alfa*(x-0.5)-0.475;
                phi2=phi1-0.05;
                theta=atand(alfa);
                if phi1<0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [1 0; 0 0.1];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    
                elseif phi1>0 && phi2<0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [100 0; 0 10];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK= [k(1,1) k(1,2) 0;k(2,1) k(2,2) 0;0 0 0];
                    
                elseif phi2>0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [1 0; 0 0.1];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0;0 0 0];
                end
                
                a=-KK*[-1;-0.2;0];
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'nikitin'
        for i=1:size(centelem,1)
            if centelem(i,2)<=(0.51-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 505 495 495 505];
            elseif (0.51-centelem(i,1))<centelem(i,2) && centelem(i,2)<(0.99-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)= [i 1000 0 0 10];
            elseif (0.99-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.49-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 10 0 0 1000];
            elseif (1.49-centelem(i,1))<=centelem(i,2)
                elem(i,5)=i;
                
                K(i,1:5)=[i 505 495 495 505];
            end
        end
    case 'nikitin1'
        for i=1:size(centelem,1)
            if centelem(i,2)<=(0.251-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 505 495 495 505];
            elseif (0.251-centelem(i,1))<centelem(i,2) && centelem(i,2)<(0.51-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)= [i 10 0 0 1000];
            elseif (0.51-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(0.751-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 1000 0 0 10];
                
                
            elseif (0.751-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(0.99-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 10 0 0 1000];
            elseif (0.99-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.249-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 1000 0 0 10];
            elseif (1.249-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.49-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 10 0 0 1000];
            elseif (1.49-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.749-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 1000 0 0 10];
            elseif (1.749-centelem(i,1))<=centelem(i,2)
                elem(i,5)=i;
                
                K(i,1:5)=[i 505 495 495 505];
                
            end
        end
    case 'lamine'
        K(1,1:5)=[1 5.5 4.5 4.5 5.5];
        elem(:,5)=1;
    case 'durlofsky'
        
        K(1,1:5)=[1 1 0 0 1];
        elem(:,5)=1;
    case 'shuec1'
        for i = 1:size(centelem,1)
            epsilon= rand(1,1);
            s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
            kmap(i,1:5)=[i s 0 0 s];
            K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
            
        end
    case 'shuec2'
        v=randperm(length(centelem));
        for i = 1:size(centelem,1)
            if centelem(i,1)<0.5
                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end
                
                kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            else
                
                epsilon= rand(1,1);
                s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                kmap(i,1:5)=[i s 0 0 s];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            end
        end
    case 'shuec3'
        v=randperm(length(centelem));
        for i = 1:size(centelem,1)
            if centelem(i,1)<0.5 && centelem(i,2)<0.5
                epsilon= rand(1,1);
                s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                kmap(i,1:5)=[i s 0 0 s];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            elseif centelem(i,1)>0.5 && centelem(i,2)<0.5
                
                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end
                
                kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            elseif centelem(i,1)<0.5 && centelem(i,2)>0.5
                
                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end
                
                kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            else
                
                epsilon= rand(1,1);
                s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                kmap(i,1:5)=[i s 0 0 s];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            end
            
        end
    case 'altamenteheterogeneo'
        % mexer m e l
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            k0=1;
            s=4;
            m=3;
            l=3;
            r=k0*exp(sqrt(s)*cos(2*pi*m*x)*cos(2*pi*l*y));
            K(i,1:5) = [i r 0 0 r];
            normKmap(i,1)= r;
            elem(i,5)=i;
            u=0;
        end
    case 'pinchout'
        % não de definimos porque tudo esta definido no start
        elem=elem;
        K=kmap;
    case 'shuecv1'
        
        K=[1	0.280963297327389	0	0	0.280963297327389
            2	0.457833366188181	0	0	0.457833366188181
            3	0.613680312903669	0	0	0.613680312903669
            4	0.376603518711135	0	0	0.376603518711135
            5	0.415237609963345	0	0	0.415237609963345
            6	0.676634555202037	0	0	0.676634555202037
            7	0.613686953306765	0	0	0.613686953306765
            8	0.376610834660443	0	0	0.376610834660443
            9	0.281020697599122	0	0	0.281020697599122
            10	0.457885474711260	0	0	0.457885474711260
            11	0.281296623870552	0	0	0.281296623870552
            12	0.172788731916259	0	0	0.172788731916259
            13	0.0889695302628119	0	0	0.0889695302628119
            14	0.143584070423681	0	0	0.143584070423681
            15	0.0664849076473011	0	0	0.0664849076473011
            16	0.0500000000000000	0	0	0.0500000000000000
            17	0.0500000000000000	0	0	0.0500000000000000
            18	0.0500000000000000	0	0	0.0500000000000000
            19	0.0848404941771196	0	0	0.0848404941771196
            20	0.0905641564150992	0	0	0.0905641564150992
            21	0.210581006378668	0	0	0.210581006378668
            22	0.192142499500473	0	0	0.192142499500473
            23	0.378602619333417	0	0	0.378602619333417
            24	0.415793864615263	0	0	0.415793864615263
            25	0.677459272624198	0	0	0.677459272624198
            26	0.617187481686624	0	0	0.617187481686624
            27	0.828234357340400	0	0	0.828234357340400
            28	0.908269893968886	0	0	0.908269893968886
            29	0.950000000000000	0	0	0.950000000000000
            30	0.914536824144816	0	0	0.914536824144816
            31	0.830930431338560	0	0	0.830930431338560
            32	0.908891105122746	0	0	0.908891105122746
            33	0.678384762047784	0	0	0.678384762047784
            34	0.621255993605384	0	0	0.621255993605384
            35	0.382255959118421	0	0	0.382255959118421
            36	0.416543820568670	0	0	0.416543820568670
            37	0.210420365544930	0	0	0.210420365544930
            38	0.193584366982287	0	0	0.193584366982287
            39	0.0807432738571369	0	0	0.0807432738571369
            40	0.0875004087123613	0	0	0.0875004087123613
            41	0.0500000000000000	0	0	0.0500000000000000
            42	0.0500000000000000	0	0	0.0500000000000000
            43	0.0500000000000000	0	0	0.0500000000000000
            44	0.0500000000000000	0	0	0.0500000000000000
            45	0.0500000000000000	0	0	0.0500000000000000
            46	0.0500000000000000	0	0	0.0500000000000000
            47	0.0500000000000000	0	0	0.0500000000000000
            48	0.0500000000000000	0	0	0.0500000000000000
            49	0.0874671909007737	0	0	0.0874671909007737
            50	0.0793444652218452	0	0	0.0793444652218452
            51	0.191937063625790	0	0	0.191937063625790
            52	0.211565459249257	0	0	0.211565459249257
            53	0.423647833486441	0	0	0.423647833486441
            54	0.384428221049007	0	0	0.384428221049007
            55	0.641370563550227	0	0	0.641370563550227
            56	0.706531155620756	0	0	0.706531155620756
            57	0.950000000000000	0	0	0.950000000000000
            58	0.903268784977887	0	0	0.903268784977887
            59	0.950000000000000	0	0	0.950000000000000
            60	0.950000000000000	0	0	0.950000000000000
            61	0.950000000000000	0	0	0.950000000000000
            62	0.950000000000000	0	0	0.950000000000000
            63	0.950000000000000	0	0	0.950000000000000
            64	0.950000000000000	0	0	0.950000000000000
            65	0.950000000000000	0	0	0.950000000000000
            66	0.950000000000000	0	0	0.950000000000000
            67	0.950000000000000	0	0	0.950000000000000
            68	0.950000000000000	0	0	0.950000000000000
            69	0.950000000000000	0	0	0.950000000000000
            70	0.907169733830307	0	0	0.907169733830307
            71	0.644111877537444	0	0	0.644111877537444
            72	0.707164658205380	0	0	0.707164658205380
            73	0.423994292374274	0	0	0.423994292374274
            74	0.385929974036454	0	0	0.385929974036454
            75	0.192587134776610	0	0	0.192587134776610
            76	0.211709864868617	0	0	0.211709864868617
            77	0.0874627557807409	0	0	0.0874627557807409
            78	0.0795215464990338	0	0	0.0795215464990338
            79	0.0500000000000000	0	0	0.0500000000000000
            80	0.0500000000000000	0	0	0.0500000000000000
            81	0.0500000000000000	0	0	0.0500000000000000
            82	0.0500000000000000	0	0	0.0500000000000000
            83	0.0500000000000000	0	0	0.0500000000000000
            84	0.0500000000000000	0	0	0.0500000000000000
            85	0.0500000000000000	0	0	0.0500000000000000
            86	0.0500000000000000	0	0	0.0500000000000000
            87	0.0500000000000000	0	0	0.0500000000000000
            88	0.0500000000000000	0	0	0.0500000000000000
            89	0.0500000000000000	0	0	0.0500000000000000
            90	0.0500000000000000	0	0	0.0500000000000000
            91	0.0500000000000000	0	0	0.0500000000000000
            92	0.0500000000000000	0	0	0.0500000000000000
            93	0.0500000000000000	0	0	0.0500000000000000
            94	0.0500000000000000	0	0	0.0500000000000000
            95	0.0500000000000000	0	0	0.0500000000000000
            96	0.0500000000000000	0	0	0.0500000000000000
            97	0.0500000000000000	0	0	0.0500000000000000
            98	0.0500000000000000	0	0	0.0500000000000000
            99	0.0500000000000000	0	0	0.0500000000000000
            100	0.0500000000000000	0	0	0.0500000000000000
            101	0.0500000000000000	0	0	0.0500000000000000
            102	0.0500000000000000	0	0	0.0500000000000000
            103	0.0500000000000000	0	0	0.0500000000000000
            104	0.0500000000000000	0	0	0.0500000000000000
            105	0.0500000000000000	0	0	0.0500000000000000
            106	0.0870383676647542	0	0	0.0870383676647542
            107	0.141830159088874	0	0	0.141830159088874
            108	0.0588931054702875	0	0	0.0588931054702875
            109	0.0789403717082000	0	0	0.0789403717082000
            110	0.190109273205522	0	0	0.190109273205522
            111	0.209611387149004	0	0	0.209611387149004
            112	0.0870383676548023	0	0	0.0870383676548023
            113	0.0789403717081389	0	0	0.0789403717081389
            114	0.190109273204901	0	0	0.190109273204901
            115	0.141830159085794	0	0	0.141830159085794
            116	0.0588931054698742	0	0	0.0588931054698742
            117	0.0500000000000000	0	0	0.0500000000000000
            118	0.0870383676550619	0	0	0.0870383676550619
            119	0.0500000000000000	0	0	0.0500000000000000
            120	0.0500000000000000	0	0	0.0500000000000000
            121	0.0500000000000000	0	0	0.0500000000000000
            122	0.0500000000000000	0	0	0.0500000000000000
            123	0.0500000000000000	0	0	0.0500000000000000
            124	0.0500000000000000	0	0	0.0500000000000000
            125	0.0500000000000000	0	0	0.0500000000000000
            126	0.0500000000000000	0	0	0.0500000000000000
            127	0.0500000000000000	0	0	0.0500000000000000
            128	0.0500000000000000	0	0	0.0500000000000000
            129	0.613680254499356	0	0	0.613680254499356
            130	0.676633848214485	0	0	0.676633848214485
            131	0.906960646567284	0	0	0.906960646567284
            132	0.822577608548544	0	0	0.822577608548544
            133	0.906961148463720	0	0	0.906961148463720
            134	0.950000000000000	0	0	0.950000000000000
            135	0.906963751456970	0	0	0.906963751456970
            136	0.822582580836316	0	0	0.822582580836316
            137	0.613719303298052	0	0	0.613719303298052
            138	0.676658319719916	0	0	0.676658319719916
            139	0.415394214915400	0	0	0.415394214915400
            140	0.376853479656303	0	0	0.376853479656303
            141	0.191426394365527	0	0	0.191426394365527
            142	0.210444925928744	0	0	0.210444925928744
            143	0.0906747893633160	0	0	0.0906747893633160
            144	0.0846495623291991	0	0	0.0846495623291991
            145	0.0500000000000000	0	0	0.0500000000000000
            146	0.0500000000000000	0	0	0.0500000000000000
            147	0.0500000000000000	0	0	0.0500000000000000
            148	0.0673526092263221	0	0	0.0673526092263221
            149	0.146168587430473	0	0	0.146168587430473
            150	0.0965444992972556	0	0	0.0965444992972556
            151	0.191032854139282	0	0	0.191032854139282
            152	0.287527726781179	0	0	0.287527726781179
            153	0.470230078580488	0	0	0.470230078580488
            154	0.317162219155413	0	0	0.317162219155413
            155	0.435503938775314	0	0	0.435503938775314
            156	0.633802740825740	0	0	0.633802740825740
            157	0.703597783315321	0	0	0.703597783315321
            158	0.494177982084070	0	0	0.494177982084070
            159	0.463641889399345	0	0	0.463641889399345
            160	0.643409532602738	0	0	0.643409532602738
            161	0.484796640772797	0	0	0.484796640772797
            162	0.359903699290704	0	0	0.359903699290704
            163	0.231315057336987	0	0	0.231315057336987
            164	0.301079615508009	0	0	0.301079615508009
            165	0.154179820437906	0	0	0.154179820437906
            166	0.123182973849001	0	0	0.123182973849001
            167	0.0544096123656189	0	0	0.0544096123656189
            168	0.0651635259759542	0	0	0.0651635259759542
            169	0.0500000000000000	0	0	0.0500000000000000
            170	0.0500000000000000	0	0	0.0500000000000000
            171	0.0500000000000000	0	0	0.0500000000000000
            172	0.0500000000000000	0	0	0.0500000000000000
            173	0.0500000000000000	0	0	0.0500000000000000
            174	0.0500000000000000	0	0	0.0500000000000000
            175	0.0500000000000000	0	0	0.0500000000000000
            176	0.0500000000000000	0	0	0.0500000000000000
            177	0.0592485624091844	0	0	0.0592485624091844
            178	0.0500000000000000	0	0	0.0500000000000000
            179	0.0885750431060262	0	0	0.0885750431060262
            180	0.143391394038956	0	0	0.143391394038956
            181	0.287499666771302	0	0	0.287499666771302
            182	0.178478575602293	0	0	0.178478575602293
            183	0.300883736531229	0	0	0.300883736531229
            184	0.480536807775710	0	0	0.480536807775710
            185	0.678804944527804	0	0	0.678804944527804
            186	0.430989231803617	0	0	0.430989231803617
            187	0.538416735932655	0	0	0.538416735932655
            188	0.830808671311961	0	0	0.830808671311961
            189	0.914759337820319	0	0	0.914759337820319
            190	0.607918182126486	0	0	0.607918182126486
            191	0.640866957647989	0	0	0.640866957647989
            192	0.942629952066226	0	0	0.942629952066226
            193	0.924372760067535	0	0	0.924372760067535
            194	0.636063442683138	0	0	0.636063442683138
            195	0.581215568191154	0	0	0.581215568191154
            196	0.845427233720139	0	0	0.845427233720139
            197	0.692689150699661	0	0	0.692689150699661
            198	0.471638091981980	0	0	0.471638091981980
            199	0.329449493761314	0	0	0.329449493761314
            200	0.490293821333525	0	0	0.490293821333525
            201	0.292846674469195	0	0	0.292846674469195
            202	0.194134254753767	0	0	0.194134254753767
            203	0.0954156407813886	0	0	0.0954156407813886
            204	0.145723952165574	0	0	0.145723952165574
            205	0.0600252545387596	0	0	0.0600252545387596
            206	0.0500000000000000	0	0	0.0500000000000000
            207	0.0500000000000000	0	0	0.0500000000000000
            208	0.0500000000000000	0	0	0.0500000000000000
            209	0.0500000000000000	0	0	0.0500000000000000
            210	0.0500000000000000	0	0	0.0500000000000000
            211	0.0500000000000000	0	0	0.0500000000000000
            212	0.0500000000000000	0	0	0.0500000000000000
            213	0.0500000000000000	0	0	0.0500000000000000
            214	0.0500000000000000	0	0	0.0500000000000000
            215	0.0500000000000000	0	0	0.0500000000000000
            216	0.0500000000000000	0	0	0.0500000000000000
            217	0.0500000000000000	0	0	0.0500000000000000
            218	0.0500000000000000	0	0	0.0500000000000000
            219	0.0500000000000000	0	0	0.0500000000000000
            220	0.0500000000000000	0	0	0.0500000000000000
            221	0.0500000000000000	0	0	0.0500000000000000
            222	0.0500000000000000	0	0	0.0500000000000000
            223	0.0500000000000000	0	0	0.0500000000000000
            224	0.0500000000000000	0	0	0.0500000000000000
            225	0.0500000000000000	0	0	0.0500000000000000
            226	0.0500000000000000	0	0	0.0500000000000000
            227	0.0500000000000000	0	0	0.0500000000000000
            228	0.0500000000000000	0	0	0.0500000000000000
            229	0.0500000000000000	0	0	0.0500000000000000
            230	0.0588931197241504	0	0	0.0588931197241504
            231	0.141830165007345	0	0	0.141830165007345
            232	0.0870383680810521	0	0	0.0870383680810521
            233	0.172421624039233	0	0	0.172421624039233
            234	0.280963294482410	0	0	0.280963294482410
            235	0.457833362340308	0	0	0.457833362340308
            236	0.280963292499216	0	0	0.280963292499216
            237	0.376603450718408	0	0	0.376603450718408
            238	0.613680251328595	0	0	0.613680251328595
            239	0.676633846184896	0	0	0.676633846184896
            240	0.415236828681174	0	0	0.415236828681174
            241	0.376603450709106	0	0	0.376603450709106
            242	0.613680251200895	0	0	0.613680251200895
            243	0.457833361770676	0	0	0.457833361770676
            244	0.280963292457247	0	0	0.280963292457247
            245	0.172421623892246	0	0	0.172421623892246
            246	0.280963292457474	0	0	0.280963292457474
            247	0.141830159085962	0	0	0.141830159085962
            248	0.0870383676551001	0	0	0.0870383676551001
            249	0.0500000000000000	0	0	0.0500000000000000
            250	0.0588931054699545	0	0	0.0588931054699545
            251	0.0500000000000000	0	0	0.0500000000000000
            252	0.0500000000000000	0	0	0.0500000000000000
            253	0.0500000000000000	0	0	0.0500000000000000
            254	0.0500000000000000	0	0	0.0500000000000000
            255	0.0500000000000000	0	0	0.0500000000000000
            256	0.0500000000000000	0	0	0.0500000000000000
            257	0.613680252422462	0	0	0.613680252422462
            258	0.457833364347552	0	0	0.457833364347552
            259	0.613680265069342	0	0	0.613680265069342
            260	0.822577578094401	0	0	0.822577578094401
            261	0.906960796897691	0	0	0.906960796897691
            262	0.676633956648424	0	0	0.676633956648424
            263	0.613681246357287	0	0	0.613681246357287
            264	0.822579275601313	0	0	0.822579275601313
            265	0.613693831718102	0	0	0.613693831718102
            266	0.457841468451352	0	0	0.457841468451352
            267	0.281019127775238	0	0	0.281019127775238
            268	0.376692411706335	0	0	0.376692411706335
            269	0.190590850040795	0	0	0.190590850040795
            270	0.142150641420937	0	0	0.142150641420937
            271	0.0604217784775940	0	0	0.0604217784775940
            272	0.0810958671683713	0	0	0.0810958671683713
            273	0.0500000000000000	0	0	0.0500000000000000
            274	0.0500000000000000	0	0	0.0500000000000000
            275	0.0500000000000000	0	0	0.0500000000000000
            276	0.0500000000000000	0	0	0.0500000000000000
            277	0.0639320135238181	0	0	0.0639320135238181
            278	0.0556919938834115	0	0	0.0556919938834115
            279	0.123428233589362	0	0	0.123428233589362
            280	0.131308153870661	0	0	0.131308153870661
            281	0.228920578731601	0	0	0.228920578731601
            282	0.231353570527768	0	0	0.231353570527768
            283	0.359908683169557	0	0	0.359908683169557
            284	0.331946129942060	0	0	0.331946129942060
            285	0.399721374817056	0	0	0.399721374817056
            286	0.463642415065552	0	0	0.463642415065552
            287	0.494177410583213	0	0	0.494177410583213
            288	0.399720737296427	0	0	0.399720737296427
            289	0.331939479196202	0	0	0.331939479196202
            290	0.435496889770342	0	0	0.435496889770342
            291	0.317105425118535	0	0	0.317105425118535
            292	0.228868736811814	0	0	0.228868736811814
            293	0.130976911504470	0	0	0.130976911504470
            294	0.190667210241458	0	0	0.190667210241458
            295	0.0946198516214560	0	0	0.0946198516214560
            296	0.0621932797796519	0	0	0.0621932797796519
            297	0.0500000000000000	0	0	0.0500000000000000
            298	0.0500000000000000	0	0	0.0500000000000000
            299	0.0500000000000000	0	0	0.0500000000000000
            300	0.0500000000000000	0	0	0.0500000000000000
            301	0.0500000000000000	0	0	0.0500000000000000
            302	0.0500000000000000	0	0	0.0500000000000000
            303	0.0500000000000000	0	0	0.0500000000000000
            304	0.0500000000000000	0	0	0.0500000000000000
            305	0.0500000000000000	0	0	0.0500000000000000
            306	0.0500000000000000	0	0	0.0500000000000000
            307	0.0500000000000000	0	0	0.0500000000000000
            308	0.0500000000000000	0	0	0.0500000000000000
            309	0.0950225692300271	0	0	0.0950225692300271
            310	0.0500000000000000	0	0	0.0500000000000000
            311	0.0976233366314859	0	0	0.0976233366314859
            312	0.166306220461475	0	0	0.166306220461475
            313	0.252290516769490	0	0	0.252290516769490
            314	0.173555154222535	0	0	0.173555154222535
            315	0.277704755171954	0	0	0.277704755171954
            316	0.340586735518262	0	0	0.340586735518262
            317	0.418977828244853	0	0	0.418977828244853
            318	0.396045524414947	0	0	0.396045524414947
            319	0.494389825396051	0	0	0.494389825396051
            320	0.473769605634859	0	0	0.473769605634859
            321	0.486759033930290	0	0	0.486759033930290
            322	0.530318804139225	0	0	0.530318804139225
            323	0.481886074568936	0	0	0.481886074568936
            324	0.443657596001147	0	0	0.443657596001147
            325	0.350183665842043	0	0	0.350183665842043
            326	0.367479348142965	0	0	0.367479348142965
            327	0.233902313394282	0	0	0.233902313394282
            328	0.235100109167083	0	0	0.235100109167083
            329	0.132726147227141	0	0	0.132726147227141
            330	0.123906919817465	0	0	0.123906919817465
            331	0.0545537950740529	0	0	0.0545537950740529
            332	0.0625858120645711	0	0	0.0625858120645711
            333	0.0500000000000000	0	0	0.0500000000000000
            334	0.0500000000000000	0	0	0.0500000000000000
            335	0.0500000000000000	0	0	0.0500000000000000
            336	0.0500000000000000	0	0	0.0500000000000000
            337	0.0500000000000000	0	0	0.0500000000000000
            338	0.0500000000000000	0	0	0.0500000000000000
            339	0.0500000000000000	0	0	0.0500000000000000
            340	0.0500000000000000	0	0	0.0500000000000000
            341	0.0500000000000000	0	0	0.0500000000000000
            342	0.0500000000000000	0	0	0.0500000000000000
            343	0.0500000000000000	0	0	0.0500000000000000
            344	0.0500000000000000	0	0	0.0500000000000000
            345	0.0500000000000000	0	0	0.0500000000000000
            346	0.0500000000000000	0	0	0.0500000000000000
            347	0.0500000000000000	0	0	0.0500000000000000
            348	0.0500000000000000	0	0	0.0500000000000000
            349	0.0500000000000000	0	0	0.0500000000000000
            350	0.0500000000000000	0	0	0.0500000000000000
            351	0.0500000000000000	0	0	0.0500000000000000
            352	0.0500000000000000	0	0	0.0500000000000000
            353	0.0500000000000000	0	0	0.0500000000000000
            354	0.0500000000000000	0	0	0.0500000000000000
            355	0.0500000000000000	0	0	0.0500000000000000
            356	0.0500000000000000	0	0	0.0500000000000000
            357	0.0789405354497030	0	0	0.0789405354497030
            358	0.0870399151082089	0	0	0.0870399151082089
            359	0.209612029713445	0	0	0.209612029713445
            360	0.190109341200208	0	0	0.190109341200208
            361	0.376603473937452	0	0	0.376603473937452
            362	0.415237048161468	0	0	0.415237048161468
            363	0.676633907829365	0	0	0.676633907829365
            364	0.613680257725459	0	0	0.613680257725459
            365	0.822577563906948	0	0	0.822577563906948
            366	0.906960632139970	0	0	0.906960632139970
            367	0.950000000000000	0	0	0.950000000000000
            368	0.906960618173464	0	0	0.906960618173464
            369	0.822577562443227	0	0	0.822577562443227
            370	0.906960618311274	0	0	0.906960618311274
            371	0.676633846215631	0	0	0.676633846215631
            372	0.613680251203314	0	0	0.613680251203314
            373	0.376603450710044	0	0	0.376603450710044
            374	0.415236828686153	0	0	0.415236828686153
            375	0.209611387149960	0	0	0.209611387149960
            376	0.190109273205310	0	0	0.190109273205310
            377	0.0789403717083057	0	0	0.0789403717083057
            378	0.0870383676550649	0	0	0.0870383676550649
            379	0.0500000000000000	0	0	0.0500000000000000
            380	0.0500000000000000	0	0	0.0500000000000000
            381	0.0500000000000000	0	0	0.0500000000000000
            382	0.0500000000000000	0	0	0.0500000000000000
            383	0.0500000000000000	0	0	0.0500000000000000
            384	0.0500000000000000	0	0	0.0500000000000000
            385	0.280963316280697	0	0	0.280963316280697
            386	0.141830383277440	0	0	0.141830383277440
            387	0.190109925613772	0	0	0.190109925613772
            388	0.376603524303269	0	0	0.376603524303269
            389	0.415237061907851	0	0	0.415237061907851
            390	0.209613007199560	0	0	0.209613007199560
            391	0.190113112097956	0	0	0.190113112097956
            392	0.376604523293273	0	0	0.376604523293273
            393	0.280970522525085	0	0	0.280970522525085
            394	0.141842331276766	0	0	0.141842331276766
            395	0.0871001809968231	0	0	0.0871001809968231
            396	0.172472075582382	0	0	0.172472075582382
            397	0.0873448624676798	0	0	0.0873448624676798
            398	0.0500000000000000	0	0	0.0500000000000000
            399	0.0500000000000000	0	0	0.0500000000000000
            400	0.0500000000000000	0	0	0.0500000000000000
            401	0.0500000000000000	0	0	0.0500000000000000
            402	0.0500000000000000	0	0	0.0500000000000000
            403	0.0500000000000000	0	0	0.0500000000000000
            404	0.0500000000000000	0	0	0.0500000000000000
            405	0.0659264498561198	0	0	0.0659264498561198
            406	0.0810974312614083	0	0	0.0810974312614083
            407	0.193655351222348	0	0	0.193655351222348
            408	0.154327297625878	0	0	0.154327297625878
            409	0.301103002182213	0	0	0.301103002182213
            410	0.382267979062808	0	0	0.382267979062808
            411	0.621259243935716	0	0	0.621259243935716
            412	0.484799881694539	0	0	0.484799881694539
            413	0.643410436019063	0	0	0.643410436019063
            414	0.830935554121262	0	0	0.830935554121262
            415	0.914548342162392	0	0	0.914548342162392
            416	0.703598650738072	0	0	0.703598650738072
            417	0.633799114837326	0	0	0.633799114837326
            418	0.828253195377148	0	0	0.828253195377148
            419	0.617187491879631	0	0	0.617187491879631
            420	0.470183046476686	0	0	0.470183046476686
            421	0.287201795817812	0	0	0.287201795817812
            422	0.378406435013078	0	0	0.378406435013078
            423	0.190894331654198	0	0	0.190894331654198
            424	0.144426890579220	0	0	0.144426890579220
            425	0.0597951301371716	0	0	0.0597951301371716
            426	0.0792449996093531	0	0	0.0792449996093531
            427	0.0500000000000000	0	0	0.0500000000000000
            428	0.0500000000000000	0	0	0.0500000000000000
            429	0.0500000000000000	0	0	0.0500000000000000
            430	0.0500000000000000	0	0	0.0500000000000000
            431	0.0500000000000000	0	0	0.0500000000000000
            432	0.0500000000000000	0	0	0.0500000000000000
            433	0.0500000000000000	0	0	0.0500000000000000
            434	0.0500000000000000	0	0	0.0500000000000000
            435	0.0500000000000000	0	0	0.0500000000000000
            436	0.0500000000000000	0	0	0.0500000000000000
            437	0.0500000000000000	0	0	0.0500000000000000
            438	0.0500000000000000	0	0	0.0500000000000000
            439	0.0848987048963925	0	0	0.0848987048963925
            440	0.0799006530691061	0	0	0.0799006530691061
            441	0.171387714456765	0	0	0.171387714456765
            442	0.198464100382300	0	0	0.198464100382300
            443	0.386745955433318	0	0	0.386745955433318
            444	0.316929090065494	0	0	0.316929090065494
            445	0.497144321476026	0	0	0.497144321476026
            446	0.624748158009121	0	0	0.624748158009121
            447	0.833893191680967	0	0	0.833893191680967
            448	0.653913386710624	0	0	0.653913386710624
            449	0.715942382711945	0	0	0.715942382711945
            450	0.918009688420598	0	0	0.918009688420598
            451	0.832682933013900	0	0	0.832682933013900
            452	0.649641269348021	0	0	0.649641269348021
            453	0.487384234412372	0	0	0.487384234412372
            454	0.621983424036578	0	0	0.621983424036578
            455	0.382503686827760	0	0	0.382503686827760
            456	0.301962927478046	0	0	0.301962927478046
            457	0.154423153397844	0	0	0.154423153397844
            458	0.193647483207236	0	0	0.193647483207236
            459	0.0807074813219627	0	0	0.0807074813219627
            460	0.0651821001927687	0	0	0.0651821001927687
            461	0.0500000000000000	0	0	0.0500000000000000
            462	0.0500000000000000	0	0	0.0500000000000000
            463	0.0500000000000000	0	0	0.0500000000000000
            464	0.0500000000000000	0	0	0.0500000000000000
            465	0.0500000000000000	0	0	0.0500000000000000
            466	0.0500000000000000	0	0	0.0500000000000000
            467	0.0500000000000000	0	0	0.0500000000000000
            468	0.0500000000000000	0	0	0.0500000000000000
            469	0.0500000000000000	0	0	0.0500000000000000
            470	0.0500000000000000	0	0	0.0500000000000000
            471	0.0500000000000000	0	0	0.0500000000000000
            472	0.0500000000000000	0	0	0.0500000000000000
            473	0.0500000000000000	0	0	0.0500000000000000
            474	0.0500000000000000	0	0	0.0500000000000000
            475	0.0500000000000000	0	0	0.0500000000000000
            476	0.0500000000000000	0	0	0.0500000000000000
            477	0.0500000000000000	0	0	0.0500000000000000
            478	0.0500000000000000	0	0	0.0500000000000000
            479	0.0500000000000000	0	0	0.0500000000000000
            480	0.0500000000000000	0	0	0.0500000000000000
            481	0.0500000000000000	0	0	0.0500000000000000
            482	0.0500000000000000	0	0	0.0500000000000000
            483	0.0500000000000000	0	0	0.0500000000000000
            484	0.0500000000000000	0	0	0.0500000000000000
            485	0.0789524016145952	0	0	0.0789524016145952
            486	0.0589700342945642	0	0	0.0589700342945642
            487	0.141862102950321	0	0	0.141862102950321
            488	0.190114268479553	0	0	0.190114268479553
            489	0.376605156925641	0	0	0.376605156925641
            490	0.280974203425115	0	0	0.280974203425115
            491	0.457836427378654	0	0	0.457836427378654
            492	0.613680730585011	0	0	0.613680730585011
            493	0.822577673192001	0	0	0.822577673192001
            494	0.613680959709574	0	0	0.613680959709574
            495	0.676633980856308	0	0	0.676633980856308
            496	0.906960638949634	0	0	0.906960638949634
            497	0.822577565691995	0	0	0.822577565691995
            498	0.613680272260744	0	0	0.613680272260744
            499	0.457833364479306	0	0	0.457833364479306
            500	0.613680251620927	0	0	0.613680251620927
            501	0.376603450753980	0	0	0.376603450753980
            502	0.280963292743705	0	0	0.280963292743705
            503	0.141830159110635	0	0	0.141830159110635
            504	0.190109273208974	0	0	0.190109273208974
            505	0.0789403717084883	0	0	0.0789403717084883
            506	0.0588931054716024	0	0	0.0588931054716024
            507	0.0500000000000000	0	0	0.0500000000000000
            508	0.0500000000000000	0	0	0.0500000000000000
            509	0.0500000000000000	0	0	0.0500000000000000
            510	0.0500000000000000	0	0	0.0500000000000000
            511	0.0500000000000000	0	0	0.0500000000000000
            512	0.0500000000000000	0	0	0.0500000000000000
            513	0.0588948736176366	0	0	0.0588948736176366
            514	0.0500000000000000	0	0	0.0500000000000000
            515	0.0500000000000000	0	0	0.0500000000000000
            516	0.0789454399082317	0	0	0.0789454399082317
            517	0.0870505286027659	0	0	0.0870505286027659
            518	0.0500000000000000	0	0	0.0500000000000000
            519	0.0500000000000000	0	0	0.0500000000000000
            520	0.0789650439252293	0	0	0.0789650439252293
            521	0.0589394595521263	0	0	0.0589394595521263
            522	0.0500000000000000	0	0	0.0500000000000000
            523	0.0500000000000000	0	0	0.0500000000000000
            524	0.0500000000000000	0	0	0.0500000000000000
            525	0.0500000000000000	0	0	0.0500000000000000
            526	0.0500000000000000	0	0	0.0500000000000000
            527	0.0500000000000000	0	0	0.0500000000000000
            528	0.0500000000000000	0	0	0.0500000000000000
            529	0.0500000000000000	0	0	0.0500000000000000
            530	0.0500000000000000	0	0	0.0500000000000000
            531	0.0500000000000000	0	0	0.0500000000000000
            532	0.0500000000000000	0	0	0.0500000000000000
            533	0.0876232008967434	0	0	0.0876232008967434
            534	0.0791523249530048	0	0	0.0791523249530048
            535	0.190308244050223	0	0	0.190308244050223
            536	0.210450629327012	0	0	0.210450629327012
            537	0.416552806143841	0	0	0.416552806143841
            538	0.376880514923595	0	0	0.376880514923595
            539	0.614073779068083	0	0	0.614073779068083
            540	0.678396724509514	0	0	0.678396724509514
            541	0.908923193669066	0	0	0.908923193669066
            542	0.823113315033596	0	0	0.823113315033596
            543	0.907698253312462	0	0	0.907698253312462
            544	0.950000000000000	0	0	0.950000000000000
            545	0.908419225832814	0	0	0.908419225832814
            546	0.823627516297096	0	0	0.823627516297096
            547	0.615138966626327	0	0	0.615138966626327
            548	0.677683787868246	0	0	0.677683787868246
            549	0.415974479738660	0	0	0.415974479738660
            550	0.378431939267073	0	0	0.378431939267073
            551	0.192075479587822	0	0	0.192075479587822
            552	0.210147295223721	0	0	0.210147295223721
            553	0.0874331657936981	0	0	0.0874331657936981
            554	0.0807190960427153	0	0	0.0807190960427153
            555	0.0500000000000000	0	0	0.0500000000000000
            556	0.0500000000000000	0	0	0.0500000000000000
            557	0.0500000000000000	0	0	0.0500000000000000
            558	0.0500000000000000	0	0	0.0500000000000000
            559	0.0500000000000000	0	0	0.0500000000000000
            560	0.0500000000000000	0	0	0.0500000000000000
            561	0.0500000000000000	0	0	0.0500000000000000
            562	0.0500000000000000	0	0	0.0500000000000000
            563	0.0500000000000000	0	0	0.0500000000000000
            564	0.0500000000000000	0	0	0.0500000000000000
            565	0.0500000000000000	0	0	0.0500000000000000
            566	0.0500000000000000	0	0	0.0500000000000000
            567	0.0811344942833787	0	0	0.0811344942833787
            568	0.0887704561353631	0	0	0.0887704561353631
            569	0.211863311326188	0	0	0.211863311326188
            570	0.192225684095960	0	0	0.192225684095960
            571	0.378353700554332	0	0	0.378353700554332
            572	0.417820276666152	0	0	0.417820276666152
            573	0.679338681454604	0	0	0.679338681454604
            574	0.614967104172258	0	0	0.614967104172258
            575	0.823478854633781	0	0	0.823478854633781
            576	0.909649966205704	0	0	0.909649966205704
            577	0.950000000000000	0	0	0.950000000000000
            578	0.907613884727292	0	0	0.907613884727292
            579	0.823078879405628	0	0	0.823078879405628
            580	0.909306627998274	0	0	0.909306627998274
            581	0.678555783154817	0	0	0.678555783154817
            582	0.614061183861678	0	0	0.614061183861678
            583	0.376866495545086	0	0	0.376866495545086
            584	0.416601154925680	0	0	0.416601154925680
            585	0.210429358936858	0	0	0.210429358936858
            586	0.190265985281344	0	0	0.190265985281344
            587	0.0790213470163169	0	0	0.0790213470163169
            588	0.0874474832700597	0	0	0.0874474832700597
            589	0.0500000000000000	0	0	0.0500000000000000
            590	0.0500000000000000	0	0	0.0500000000000000
            591	0.0500000000000000	0	0	0.0500000000000000
            592	0.0500000000000000	0	0	0.0500000000000000
            593	0.0500000000000000	0	0	0.0500000000000000
            594	0.0500000000000000	0	0	0.0500000000000000
            595	0.0500000000000000	0	0	0.0500000000000000
            596	0.0500000000000000	0	0	0.0500000000000000
            597	0.0500000000000000	0	0	0.0500000000000000
            598	0.0500000000000000	0	0	0.0500000000000000
            599	0.0500000000000000	0	0	0.0500000000000000
            600	0.0500000000000000	0	0	0.0500000000000000
            601	0.0500000000000000	0	0	0.0500000000000000
            602	0.0500000000000000	0	0	0.0500000000000000
            603	0.0500000000000000	0	0	0.0500000000000000
            604	0.0500000000000000	0	0	0.0500000000000000
            605	0.0500000000000000	0	0	0.0500000000000000
            606	0.0500000000000000	0	0	0.0500000000000000
            607	0.0500000000000000	0	0	0.0500000000000000
            608	0.0500000000000000	0	0	0.0500000000000000
            609	0.0500000000000000	0	0	0.0500000000000000
            610	0.0500000000000000	0	0	0.0500000000000000
            611	0.0500000000000000	0	0	0.0500000000000000
            612	0.0500000000000000	0	0	0.0500000000000000
            613	0.0500000000000000	0	0	0.0500000000000000
            614	0.0500000000000000	0	0	0.0500000000000000
            615	0.0500000000000000	0	0	0.0500000000000000
            616	0.0872064054895917	0	0	0.0872064054895917
            617	0.172479020701617	0	0	0.172479020701617
            618	0.0872867584011828	0	0	0.0872867584011828
            619	0.141899951887418	0	0	0.141899951887418
            620	0.280979419213202	0	0	0.280979419213202
            621	0.376607177930335	0	0	0.376607177930335
            622	0.190125404463061	0	0	0.190125404463061
            623	0.209614454118174	0	0	0.209614454118174
            624	0.415237537281968	0	0	0.415237537281968
            625	0.376603561523954	0	0	0.376603561523954
            626	0.190109752867051	0	0	0.190109752867051
            627	0.141830220794523	0	0	0.141830220794523
            628	0.280963306712135	0	0	0.280963306712135
            629	0.172421625400320	0	0	0.172421625400320
            630	0.0870383741853611	0	0	0.0870383741853611
            631	0.0500000000000000	0	0	0.0500000000000000
            632	0.0870383677861482	0	0	0.0870383677861482
            633	0.0500000000000000	0	0	0.0500000000000000
            634	0.0500000000000000	0	0	0.0500000000000000
            635	0.0500000000000000	0	0	0.0500000000000000
            636	0.0500000000000000	0	0	0.0500000000000000
            637	0.0500000000000000	0	0	0.0500000000000000
            638	0.0500000000000000	0	0	0.0500000000000000
            639	0.0500000000000000	0	0	0.0500000000000000
            640	0.0500000000000000	0	0	0.0500000000000000
            641	0.0500000000000000	0	0	0.0500000000000000
            642	0.0500000000000000	0	0	0.0500000000000000
            643	0.0500000000000000	0	0	0.0500000000000000
            644	0.0500000000000000	0	0	0.0500000000000000
            645	0.0500000000000000	0	0	0.0500000000000000
            646	0.0500000000000000	0	0	0.0500000000000000
            647	0.0500000000000000	0	0	0.0500000000000000
            648	0.0500000000000000	0	0	0.0500000000000000
            649	0.0500000000000000	0	0	0.0500000000000000
            650	0.0500000000000000	0	0	0.0500000000000000
            651	0.0500000000000000	0	0	0.0500000000000000
            652	0.0500000000000000	0	0	0.0500000000000000
            653	0.0500000000000000	0	0	0.0500000000000000
            654	0.0500000000000000	0	0	0.0500000000000000
            655	0.0500000000000000	0	0	0.0500000000000000
            656	0.0500000000000000	0	0	0.0500000000000000
            657	0.0500000000000000	0	0	0.0500000000000000
            658	0.0500000000000000	0	0	0.0500000000000000
            659	0.0500000000000000	0	0	0.0500000000000000
            660	0.0500000000000000	0	0	0.0500000000000000
            661	0.0593289733268074	0	0	0.0593289733268074
            662	0.0500000000000000	0	0	0.0500000000000000
            663	0.0878265514601492	0	0	0.0878265514601492
            664	0.142041204091649	0	0	0.142041204091649
            665	0.281129764973877	0	0	0.281129764973877
            666	0.172923375823372	0	0	0.172923375823372
            667	0.281923598366943	0	0	0.281923598366943
            668	0.458149922239503	0	0	0.458149922239503
            669	0.614468391522712	0	0	0.614468391522712
            670	0.379214543884258	0	0	0.379214543884258
            671	0.421478195367360	0	0	0.421478195367360
            672	0.678437458879176	0	0	0.678437458879176
            673	0.617187654197123	0	0	0.617187654197123
            674	0.388953810317783	0	0	0.388953810317783
            675	0.301083919922915	0	0	0.301083919922915
            676	0.463509486170358	0	0	0.463509486170358
            677	0.288554070372854	0	0	0.288554070372854
            678	0.199397318940837	0	0	0.199397318940837
            679	0.116825646709922	0	0	0.116825646709922
            680	0.150204146976274	0	0	0.150204146976274
            681	0.0665408083051782	0	0	0.0665408083051782
            682	0.0633536335574821	0	0	0.0633536335574821
            683	0.0500000000000000	0	0	0.0500000000000000
            684	0.0500000000000000	0	0	0.0500000000000000
            685	0.0500000000000000	0	0	0.0500000000000000
            686	0.0500000000000000	0	0	0.0500000000000000
            687	0.0500000000000000	0	0	0.0500000000000000
            688	0.0500000000000000	0	0	0.0500000000000000
            689	0.0500000000000000	0	0	0.0500000000000000
            690	0.0500000000000000	0	0	0.0500000000000000
            691	0.0500000000000000	0	0	0.0500000000000000
            692	0.0500000000000000	0	0	0.0500000000000000
            693	0.0500000000000000	0	0	0.0500000000000000
            694	0.0500000000000000	0	0	0.0500000000000000
            695	0.0659336820795137	0	0	0.0659336820795137
            696	0.0673031160808080	0	0	0.0673031160808080
            697	0.149465994481747	0	0	0.149465994481747
            698	0.114019873032639	0	0	0.114019873032639
            699	0.192548056055798	0	0	0.192548056055798
            700	0.286684995525289	0	0	0.286684995525289
            701	0.461377703641897	0	0	0.461377703641897
            702	0.293317939051451	0	0	0.293317939051451
            703	0.382845049961463	0	0	0.382845049961463
            704	0.615508769275741	0	0	0.615508769275741
            705	0.677436747359362	0	0	0.677436747359362
            706	0.417834195699044	0	0	0.417834195699044
            707	0.377496323755779	0	0	0.377496323755779
            708	0.613998010934320	0	0	0.613998010934320
            709	0.457960262726404	0	0	0.457960262726404
            710	0.281219500726244	0	0	0.281219500726244
            711	0.172486578254887	0	0	0.172486578254887
            712	0.281020394953351	0	0	0.281020394953351
            713	0.141859883533906	0	0	0.141859883533906
            714	0.0870623961657947	0	0	0.0870623961657947
            715	0.0500000000000000	0	0	0.0500000000000000
            716	0.0589181838300314	0	0	0.0589181838300314
            717	0.0500000000000000	0	0	0.0500000000000000
            718	0.0500000000000000	0	0	0.0500000000000000
            719	0.0500000000000000	0	0	0.0500000000000000
            720	0.0500000000000000	0	0	0.0500000000000000
            721	0.0500000000000000	0	0	0.0500000000000000
            722	0.0500000000000000	0	0	0.0500000000000000
            723	0.0500000000000000	0	0	0.0500000000000000
            724	0.0500000000000000	0	0	0.0500000000000000
            725	0.0500000000000000	0	0	0.0500000000000000
            726	0.0500000000000000	0	0	0.0500000000000000
            727	0.0500000000000000	0	0	0.0500000000000000
            728	0.0500000000000000	0	0	0.0500000000000000
            729	0.0500000000000000	0	0	0.0500000000000000
            730	0.0589413109822564	0	0	0.0589413109822564
            731	0.0789929685340508	0	0	0.0789929685340508
            732	0.0500000000000000	0	0	0.0500000000000000
            733	0.0500000000000000	0	0	0.0500000000000000
            734	0.0871063267024143	0	0	0.0871063267024143
            735	0.0790568597913653	0	0	0.0790568597913653
            736	0.0500000000000000	0	0	0.0500000000000000
            737	0.0500000000000000	0	0	0.0500000000000000
            738	0.0591939588770976	0	0	0.0591939588770976
            739	0.0500000000000000	0	0	0.0500000000000000
            740	0.0500000000000000	0	0	0.0500000000000000
            741	0.0500000000000000	0	0	0.0500000000000000
            742	0.0500000000000000	0	0	0.0500000000000000
            743	0.0500000000000000	0	0	0.0500000000000000
            744	0.0500000000000000	0	0	0.0500000000000000
            745	0.0500000000000000	0	0	0.0500000000000000
            746	0.0500000000000000	0	0	0.0500000000000000
            747	0.0500000000000000	0	0	0.0500000000000000
            748	0.0591416476994054	0	0	0.0591416476994054
            749	0.0789978242617288	0	0	0.0789978242617288
            750	0.0500000000000000	0	0	0.0500000000000000
            751	0.0500000000000000	0	0	0.0500000000000000
            752	0.0870492924464703	0	0	0.0870492924464703
            753	0.0789420806035381	0	0	0.0789420806035381
            754	0.0500000000000000	0	0	0.0500000000000000
            755	0.0500000000000000	0	0	0.0500000000000000
            756	0.0588933253675941	0	0	0.0588933253675941
            757	0.0500000000000000	0	0	0.0500000000000000
            758	0.0500000000000000	0	0	0.0500000000000000
            759	0.0500000000000000	0	0	0.0500000000000000
            760	0.0500000000000000	0	0	0.0500000000000000
            761	0.0500000000000000	0	0	0.0500000000000000
            762	0.0500000000000000	0	0	0.0500000000000000
            763	0.0500000000000000	0	0	0.0500000000000000
            764	0.0500000000000000	0	0	0.0500000000000000
            765	0.0500000000000000	0	0	0.0500000000000000
            766	0.0500000000000000	0	0	0.0500000000000000
            767	0.0500000000000000	0	0	0.0500000000000000
            768	0.0500000000000000	0	0	0.0500000000000000
            769	0.0500000000000000	0	0	0.0500000000000000
            770	0.0500000000000000	0	0	0.0500000000000000
            771	0.0500000000000000	0	0	0.0500000000000000
            772	0.0500000000000000	0	0	0.0500000000000000
            773	0.0500000000000000	0	0	0.0500000000000000
            774	0.0500000000000000	0	0	0.0500000000000000
            775	0.0500000000000000	0	0	0.0500000000000000
            776	0.0500000000000000	0	0	0.0500000000000000
            777	0.0500000000000000	0	0	0.0500000000000000
            778	0.0597341763621223	0	0	0.0597341763621223
            779	0.0793738439446904	0	0	0.0793738439446904
            780	0.0500000000000000	0	0	0.0500000000000000
            781	0.0500000000000000	0	0	0.0500000000000000
            782	0.0872503782406056	0	0	0.0872503782406056
            783	0.0791709456907427	0	0	0.0791709456907427
            784	0.0500000000000000	0	0	0.0500000000000000
            785	0.0500000000000000	0	0	0.0500000000000000
            786	0.0596395015447239	0	0	0.0596395015447239
            787	0.0500000000000000	0	0	0.0500000000000000
            788	0.0500000000000000	0	0	0.0500000000000000
            789	0.0500000000000000	0	0	0.0500000000000000
            790	0.0500000000000000	0	0	0.0500000000000000
            791	0.0500000000000000	0	0	0.0500000000000000
            792	0.0500000000000000	0	0	0.0500000000000000
            793	0.0886497620005134	0	0	0.0886497620005134
            794	0.0500000000000000	0	0	0.0500000000000000
            795	0.0658517680620368	0	0	0.0658517680620368
            796	0.144666776193590	0	0	0.144666776193590
            797	0.197743128918615	0	0	0.197743128918615
            798	0.0973527106841748	0	0	0.0973527106841748
            799	0.131007469536497	0	0	0.131007469536497
            800	0.227867307159929	0	0	0.227867307159929
            801	0.226253686265216	0	0	0.226253686265216
            802	0.165985328028155	0	0	0.165985328028155
            803	0.200735965160787	0	0	0.200735965160787
            804	0.200728803664982	0	0	0.200728803664982
            805	0.166010870193122	0	0	0.166010870193122
            806	0.226327818440590	0	0	0.226327818440590
            807	0.228260283070481	0	0	0.228260283070481
            808	0.131143409164648	0	0	0.131143409164648
            809	0.0979118077710187	0	0	0.0979118077710187
            810	0.199435867001400	0	0	0.199435867001400
            811	0.150649967131426	0	0	0.150649967131426
            812	0.0677138257083818	0	0	0.0677138257083818
            813	0.0500000000000000	0	0	0.0500000000000000
            814	0.106014612244712	0	0	0.106014612244712
            815	0.0880738357197239	0	0	0.0880738357197239
            816	0.0500000000000000	0	0	0.0500000000000000
            817	0.0500000000000000	0	0	0.0500000000000000
            818	0.105482573002567	0	0	0.105482573002567
            819	0.148793766656048	0	0	0.148793766656048
            820	0.0632437288482464	0	0	0.0632437288482464
            821	0.0859008480150019	0	0	0.0859008480150019
            822	0.194448372905546	0	0	0.194448372905546
            823	0.217591825789563	0	0	0.217591825789563
            824	0.105451339726774	0	0	0.105451339726774
            825	0.122910038909744	0	0	0.122910038909744
            826	0.208430506615920	0	0	0.208430506615920
            827	0.177983808802560	0	0	0.177983808802560
            828	0.145937414529013	0	0	0.145937414529013
            829	0.177973373006427	0	0	0.177973373006427
            830	0.145933124747844	0	0	0.145933124747844
            831	0.122877671988878	0	0	0.122877671988878
            832	0.208354652404486	0	0	0.208354652404486
            833	0.217188239053173	0	0	0.217188239053173
            834	0.105283415276234	0	0	0.105283415276234
            835	0.0851757921511244	0	0	0.0851757921511244
            836	0.192698523707692	0	0	0.192698523707692
            837	0.142561180487314	0	0	0.142561180487314
            838	0.0606604682154218	0	0	0.0606604682154218
            839	0.0500000000000000	0	0	0.0500000000000000
            840	0.0872194012804870	0	0	0.0872194012804870
            841	0.0500000000000000	0	0	0.0500000000000000
            842	0.0500000000000000	0	0	0.0500000000000000
            843	0.0500000000000000	0	0	0.0500000000000000
            844	0.0500000000000000	0	0	0.0500000000000000
            845	0.0500000000000000	0	0	0.0500000000000000
            846	0.0500000000000000	0	0	0.0500000000000000
            847	0.0500000000000000	0	0	0.0500000000000000
            848	0.0500000000000000	0	0	0.0500000000000000
            849	0.0500000000000000	0	0	0.0500000000000000
            850	0.0500000000000000	0	0	0.0500000000000000
            851	0.0500000000000000	0	0	0.0500000000000000
            852	0.0500000000000000	0	0	0.0500000000000000
            853	0.0500000000000000	0	0	0.0500000000000000
            854	0.0889622325358130	0	0	0.0889622325358130
            855	0.174133123741361	0	0	0.174133123741361
            856	0.0873651712190774	0	0	0.0873651712190774
            857	0.142135166028813	0	0	0.142135166028813
            858	0.282566566008368	0	0	0.282566566008368
            859	0.378314423196738	0	0	0.378314423196738
            860	0.190435016381720	0	0	0.190435016381720
            861	0.209979248508818	0	0	0.209979248508818
            862	0.417156407036760	0	0	0.417156407036760
            863	0.378595022437788	0	0	0.378595022437788
            864	0.190503305772951	0	0	0.190503305772951
            865	0.142235846168232	0	0	0.142235846168232
            866	0.282746387862579	0	0	0.282746387862579
            867	0.173788952162610	0	0	0.173788952162610
            868	0.0875357216319172	0	0	0.0875357216319172
            869	0.0500000000000000	0	0	0.0500000000000000
            870	0.0880086955528293	0	0	0.0880086955528293
            871	0.0500000000000000	0	0	0.0500000000000000
            872	0.0500000000000000	0	0	0.0500000000000000
            873	0.0500000000000000	0	0	0.0500000000000000
            874	0.0500000000000000	0	0	0.0500000000000000
            875	0.0500000000000000	0	0	0.0500000000000000
            876	0.0500000000000000	0	0	0.0500000000000000
            877	0.0500000000000000	0	0	0.0500000000000000
            878	0.0500000000000000	0	0	0.0500000000000000
            879	0.0500000000000000	0	0	0.0500000000000000
            880	0.0500000000000000	0	0	0.0500000000000000
            881	0.0500000000000000	0	0	0.0500000000000000
            882	0.0500000000000000	0	0	0.0500000000000000
            883	0.0500000000000000	0	0	0.0500000000000000
            884	0.0500000000000000	0	0	0.0500000000000000
            885	0.0500000000000000	0	0	0.0500000000000000
            886	0.0500000000000000	0	0	0.0500000000000000
            887	0.0500000000000000	0	0	0.0500000000000000
            888	0.0500000000000000	0	0	0.0500000000000000
            889	0.0500000000000000	0	0	0.0500000000000000
            890	0.0500000000000000	0	0	0.0500000000000000
            891	0.0500000000000000	0	0	0.0500000000000000
            892	0.0500000000000000	0	0	0.0500000000000000
            893	0.0500000000000000	0	0	0.0500000000000000
            894	0.0500000000000000	0	0	0.0500000000000000
            895	0.0500000000000000	0	0	0.0500000000000000
            896	0.0500000000000000	0	0	0.0500000000000000
            897	0.0500000000000000	0	0	0.0500000000000000
            898	0.0500000000000000	0	0	0.0500000000000000
            899	0.0658714606580088	0	0	0.0658714606580088
            900	0.0500000000000000	0	0	0.0500000000000000
            901	0.0515199955299591	0	0	0.0515199955299591
            902	0.114002377836854	0	0	0.114002377836854
            903	0.192538132926481	0	0	0.192538132926481
            904	0.0926968861883383	0	0	0.0926968861883383
            905	0.145303771762696	0	0	0.145303771762696
            906	0.293308548401159	0	0	0.293308548401159
            907	0.382835840029187	0	0	0.382835840029187
            908	0.191864892754620	0	0	0.191864892754620
            909	0.210350889066812	0	0	0.210350889066812
            910	0.417827686353387	0	0	0.417827686353387
            911	0.377503998296541	0	0	0.377503998296541
            912	0.190415953804087	0	0	0.190415953804087
            913	0.142151536494733	0	0	0.142151536494733
            914	0.281286403226940	0	0	0.281286403226940
            915	0.172766176699402	0	0	0.172766176699402
            916	0.0879622115950832	0	0	0.0879622115950832
            917	0.0500000000000000	0	0	0.0500000000000000
            918	0.0880246182032136	0	0	0.0880246182032136
            919	0.0500000000000000	0	0	0.0500000000000000
            920	0.0500000000000000	0	0	0.0500000000000000
            921	0.0500000000000000	0	0	0.0500000000000000
            922	0.0500000000000000	0	0	0.0500000000000000
            923	0.0500000000000000	0	0	0.0500000000000000
            924	0.0500000000000000	0	0	0.0500000000000000
            925	0.0635094430732130	0	0	0.0635094430732130
            926	0.0672704270257035	0	0	0.0672704270257035
            927	0.150335897752525	0	0	0.150335897752525
            928	0.116844839684072	0	0	0.116844839684072
            929	0.199399950958334	0	0	0.199399950958334
            930	0.288567797147707	0	0	0.288567797147707
            931	0.463527061722124	0	0	0.463527061722124
            932	0.301104466889527	0	0	0.301104466889527
            933	0.389100672107832	0	0	0.389100672107832
            934	0.617397299600958	0	0	0.617397299600958
            935	0.679690935945121	0	0	0.679690935945121
            936	0.422270023581075	0	0	0.422270023581075
            937	0.382659455731192	0	0	0.382659455731192
            938	0.620059080149079	0	0	0.620059080149079
            939	0.478197544941516	0	0	0.478197544941516
            940	0.294192004150509	0	0	0.294192004150509
            941	0.208813188846828	0	0	0.208813188846828
            942	0.339926650443429	0	0	0.339926650443429
            943	0.283679510683961	0	0	0.283679510683961
            944	0.174145037020042	0	0	0.174145037020042
            945	0.208631461579458	0	0	0.208631461579458
            946	0.339875589428131	0	0	0.339875589428131
            947	0.478019405429148	0	0	0.478019405429148
            948	0.293557987989541	0	0	0.293557987989541
            949	0.380955896442564	0	0	0.380955896442564
            950	0.619580435725419	0	0	0.619580435725419
            951	0.678667079430544	0	0	0.678667079430544
            952	0.418626033508399	0	0	0.418626033508399
            953	0.382987431192427	0	0	0.382987431192427
            954	0.615679422731518	0	0	0.615679422731518
            955	0.461340615881141	0	0	0.461340615881141
            956	0.293331801358298	0	0	0.293331801358298
            957	0.192540551498789	0	0	0.192540551498789
            958	0.286620307257654	0	0	0.286620307257654
            959	0.149408071917260	0	0	0.149408071917260
            960	0.114002406637429	0	0	0.114002406637429
            961	0.0658738555702941	0	0	0.0658738555702941
            962	0.0672568853352944	0	0	0.0672568853352944
            963	0.0500000000000000	0	0	0.0500000000000000
            964	0.0500000000000000	0	0	0.0500000000000000
            965	0.0500000000000000	0	0	0.0500000000000000
            966	0.0500000000000000	0	0	0.0500000000000000
            967	0.0500000000000000	0	0	0.0500000000000000
            968	0.0500000000000000	0	0	0.0500000000000000
            969	0.0500000000000000	0	0	0.0500000000000000
            970	0.0500000000000000	0	0	0.0500000000000000
            971	0.0500000000000000	0	0	0.0500000000000000
            972	0.0500000000000000	0	0	0.0500000000000000
            973	0.0500000000000000	0	0	0.0500000000000000
            974	0.0500000000000000	0	0	0.0500000000000000
            975	0.0500000000000000	0	0	0.0500000000000000
            976	0.0500000000000000	0	0	0.0500000000000000
            977	0.0500000000000000	0	0	0.0500000000000000
            978	0.0541909504211117	0	0	0.0541909504211117
            979	0.109556483169600	0	0	0.109556483169600
            980	0.0675051593122510	0	0	0.0675051593122510
            981	0.150134558555149	0	0	0.150134558555149
            982	0.219660563728113	0	0	0.219660563728113
            983	0.402950888783117	0	0	0.402950888783117
            984	0.288366157826308	0	0	0.288366157826308
            985	0.464770183458803	0	0	0.464770183458803
            986	0.638369589355510	0	0	0.638369589355510
            987	0.848924943150015	0	0	0.848924943150015
            988	0.621082919918694	0	0	0.621082919918694
            989	0.684936674270421	0	0	0.684936674270421
            990	0.936511451945436	0	0	0.936511451945436
            991	0.853190721805242	0	0	0.853190721805242
            992	0.622282101368624	0	0	0.622282101368624
            993	0.465481875275996	0	0	0.465481875275996
            994	0.640892296040516	0	0	0.640892296040516
            995	0.396778370047588	0	0	0.396778370047588
            996	0.286642129475543	0	0	0.286642129475543
            997	0.145333568556298	0	0	0.145333568556298
            998	0.202469890494630	0	0	0.202469890494630
            999	0.0851857516309516	0	0	0.0851857516309516
            1000	0.0607213619379723	0	0	0.0607213619379723
            1001	0.0500000000000000	0	0	0.0500000000000000
            1002	0.0500000000000000	0	0	0.0500000000000000
            1003	0.0500000000000000	0	0	0.0500000000000000
            1004	0.0500000000000000	0	0	0.0500000000000000
            1005	0.0500000000000000	0	0	0.0500000000000000
            1006	0.0500000000000000	0	0	0.0500000000000000
            1007	0.0500000000000000	0	0	0.0500000000000000
            1008	0.0500000000000000	0	0	0.0500000000000000
            1009	0.0500000000000000	0	0	0.0500000000000000
            1010	0.0500000000000000	0	0	0.0500000000000000
            1011	0.0500000000000000	0	0	0.0500000000000000
            1012	0.0500000000000000	0	0	0.0500000000000000
            1013	0.0500000000000000	0	0	0.0500000000000000
            1014	0.0500000000000000	0	0	0.0500000000000000
            1015	0.0500000000000000	0	0	0.0500000000000000
            1016	0.0500000000000000	0	0	0.0500000000000000
            1017	0.0500000000000000	0	0	0.0500000000000000
            1018	0.0500000000000000	0	0	0.0500000000000000
            1019	0.0500000000000000	0	0	0.0500000000000000
            1020	0.0500000000000000	0	0	0.0500000000000000
            1021	0.0500000000000000	0	0	0.0500000000000000
            1022	0.0500000000000000	0	0	0.0500000000000000
            1023	0.0500000000000000	0	0	0.0500000000000000
            1024	0.0500000000000000	0	0	0.0500000000000000
            1025	0.0990562117508940	0	0	0.0990562117508940
            1026	0.217072505158281	0	0	0.217072505158281
            1027	0.288551763537979	0	0	0.288551763537979
            1028	0.145931534813686	0	0	0.145931534813686
            1029	0.220770598932295	0	0	0.220770598932295
            1030	0.380218551593754	0	0	0.380218551593754
            1031	0.518433614773225	0	0	0.518433614773225
            1032	0.339856461126894	0	0	0.339856461126894
            1033	0.493974958055927	0	0	0.493974958055927
            1034	0.700718625337749	0	0	0.700718625337749
            1035	0.866514524092504	0	0	0.866514524092504
            1036	0.631924627621904	0	0	0.631924627621904
            1037	0.684210307997130	0	0	0.684210307997130
            1038	0.925205054624952	0	0	0.925205054624952
            1039	0.828810634753379	0	0	0.828810634753379
            1040	0.616272274708887	0	0	0.616272274708887
            1041	0.458583199799264	0	0	0.458583199799264
            1042	0.615441545129806	0	0	0.615441545129806
            1043	0.377076397356808	0	0	0.377076397356808
            1044	0.281240758769847	0	0	0.281240758769847
            1045	0.142362977679262	0	0	0.142362977679262
            1046	0.190585140133007	0	0	0.190585140133007
            1047	0.0809597827303869	0	0	0.0809597827303869
            1048	0.0610535917934163	0	0	0.0610535917934163
            1049	0.0500000000000000	0	0	0.0500000000000000
            1050	0.0500000000000000	0	0	0.0500000000000000
            1051	0.0500000000000000	0	0	0.0500000000000000
            1052	0.0500000000000000	0	0	0.0500000000000000
            1053	0.0819974978106939	0	0	0.0819974978106939
            1054	0.0891221710583191	0	0	0.0891221710583191
            1055	0.210312086710441	0	0	0.210312086710441
            1056	0.192288764925630	0	0	0.192288764925630
            1057	0.378399839112982	0	0	0.378399839112982
            1058	0.415629278800004	0	0	0.415629278800004
            1059	0.676946581767824	0	0	0.676946581767824
            1060	0.615043536673812	0	0	0.615043536673812
            1061	0.823712585081065	0	0	0.823712585081065
            1062	0.907481248829754	0	0	0.907481248829754
            1063	0.950000000000000	0	0	0.950000000000000
            1064	0.909116217754278	0	0	0.909116217754278
            1065	0.830321338232634	0	0	0.830321338232634
            1066	0.915346193327428	0	0	0.915346193327428
            1067	0.706374730531353	0	0	0.706374730531353
            1068	0.640700949951967	0	0	0.640700949951967
            1069	0.455560108370164	0	0	0.455560108370164
            1070	0.502278873087674	0	0	0.502278873087674
            1071	0.419224009809832	0	0	0.419224009809832
            1072	0.380223014635968	0	0	0.380223014635968
            1073	0.455548276210207	0	0	0.455548276210207
            1074	0.502276264895867	0	0	0.502276264895867
            1075	0.706366263067067	0	0	0.706366263067067
            1076	0.640659724399092	0	0	0.640659724399092
            1077	0.830210649056164	0	0	0.830210649056164
            1078	0.915324401169568	0	0	0.915324401169568
            1079	0.950000000000000	0	0	0.950000000000000
            1080	0.908879490870074	0	0	0.908879490870074
            1081	0.823315066222144	0	0	0.823315066222144
            1082	0.907404563621507	0	0	0.907404563621507
            1083	0.676843793865013	0	0	0.676843793865013
            1084	0.614533966312621	0	0	0.614533966312621
            1085	0.377916999532146	0	0	0.377916999532146
            1086	0.415494081940099	0	0	0.415494081940099
            1087	0.209957050520193	0	0	0.209957050520193
            1088	0.191865823078701	0	0	0.191865823078701
            1089	0.0809028352621475	0	0	0.0809028352621475
            1090	0.0874823641787062	0	0	0.0874823641787062
            1091	0.0500000000000000	0	0	0.0500000000000000
            1092	0.0500000000000000	0	0	0.0500000000000000
            1093	0.0500000000000000	0	0	0.0500000000000000
            1094	0.0500000000000000	0	0	0.0500000000000000
            1095	0.0500000000000000	0	0	0.0500000000000000
            1096	0.0500000000000000	0	0	0.0500000000000000
            1097	0.0500000000000000	0	0	0.0500000000000000
            1098	0.0500000000000000	0	0	0.0500000000000000
            1099	0.0500000000000000	0	0	0.0500000000000000
            1100	0.0500000000000000	0	0	0.0500000000000000
            1101	0.0500000000000000	0	0	0.0500000000000000
            1102	0.0888844337364709	0	0	0.0888844337364709
            1103	0.149814329332126	0	0	0.149814329332126
            1104	0.0674307808308675	0	0	0.0674307808308675
            1105	0.109400392222228	0	0	0.109400392222228
            1106	0.218824126034270	0	0	0.218824126034270
            1107	0.294783484049794	0	0	0.294783484049794
            1108	0.176665030865657	0	0	0.176665030865657
            1109	0.296127547662969	0	0	0.296127547662969
            1110	0.398462805582887	0	0	0.398462805582887
            1111	0.562370553725281	0	0	0.562370553725281
            1112	0.492374181978637	0	0	0.492374181978637
            1113	0.748916921013239	0	0	0.748916921013239
            1114	0.787757005472150	0	0	0.787757005472150
            1115	0.950000000000000	0	0	0.950000000000000
            1116	0.950000000000000	0	0	0.950000000000000
            1117	0.950000000000000	0	0	0.950000000000000
            1118	0.950000000000000	0	0	0.950000000000000
            1119	0.950000000000000	0	0	0.950000000000000
            1120	0.950000000000000	0	0	0.950000000000000
            1121	0.756301458988751	0	0	0.756301458988751
            1122	0.805541924674090	0	0	0.805541924674090
            1123	0.518841343798702	0	0	0.518841343798702
            1124	0.474298501576948	0	0	0.474298501576948
            1125	0.245786036863057	0	0	0.245786036863057
            1126	0.277229630904606	0	0	0.277229630904606
            1127	0.122896194449242	0	0	0.122896194449242
            1128	0.105289851762825	0	0	0.105289851762825
            1129	0.0500000000000000	0	0	0.0500000000000000
            1130	0.0500000000000000	0	0	0.0500000000000000
            1131	0.0500000000000000	0	0	0.0500000000000000
            1132	0.0500000000000000	0	0	0.0500000000000000
            1133	0.0500000000000000	0	0	0.0500000000000000
            1134	0.0500000000000000	0	0	0.0500000000000000
            1135	0.0500000000000000	0	0	0.0500000000000000
            1136	0.0500000000000000	0	0	0.0500000000000000
            1137	0.0500000000000000	0	0	0.0500000000000000
            1138	0.0500000000000000	0	0	0.0500000000000000
            1139	0.0500000000000000	0	0	0.0500000000000000
            1140	0.0500000000000000	0	0	0.0500000000000000
            1141	0.0500000000000000	0	0	0.0500000000000000
            1142	0.0500000000000000	0	0	0.0500000000000000
            1143	0.0500000000000000	0	0	0.0500000000000000
            1144	0.0500000000000000	0	0	0.0500000000000000
            1145	0.0500000000000000	0	0	0.0500000000000000
            1146	0.0500000000000000	0	0	0.0500000000000000
            1147	0.0500000000000000	0	0	0.0500000000000000
            1148	0.0500000000000000	0	0	0.0500000000000000
            1149	0.0500000000000000	0	0	0.0500000000000000
            1150	0.0500000000000000	0	0	0.0500000000000000
            1151	0.0500000000000000	0	0	0.0500000000000000
            1152	0.0500000000000000	0	0	0.0500000000000000
            1153	0.406332670585499	0	0	0.406332670585499
            1154	0.640643525675154	0	0	0.640643525675154
            1155	0.755574279548658	0	0	0.755574279548658
            1156	0.502275201014243	0	0	0.502275201014243
            1157	0.586214843021296	0	0	0.586214843021296
            1158	0.803789592410240	0	0	0.803789592410240
            1159	0.834436874174203	0	0	0.834436874174203
            1160	0.696200125872794	0	0	0.696200125872794
            1161	0.849055474226208	0	0	0.849055474226208
            1162	0.894643589919108	0	0	0.894643589919108
            1163	0.950000000000000	0	0	0.950000000000000
            1164	0.950000000000000	0	0	0.950000000000000
            1165	0.950000000000000	0	0	0.950000000000000
            1166	0.950000000000000	0	0	0.950000000000000
            1167	0.842693905151482	0	0	0.842693905151482
            1168	0.919306103129214	0	0	0.919306103129214
            1169	0.680109419150155	0	0	0.680109419150155
            1170	0.619337109857933	0	0	0.619337109857933
            1171	0.377948714401072	0	0	0.377948714401072
            1172	0.416092181893260	0	0	0.416092181893260
            1173	0.210101625597241	0	0	0.210101625597241
            1174	0.190606564687080	0	0	0.190606564687080
            1175	0.0802869834830069	0	0	0.0802869834830069
            1176	0.0888250768504359	0	0	0.0888250768504359
            1177	0.0500000000000000	0	0	0.0500000000000000
            1178	0.0500000000000000	0	0	0.0500000000000000
            1179	0.0500000000000000	0	0	0.0500000000000000
            1180	0.0500000000000000	0	0	0.0500000000000000
            1181	0.0809231867781842	0	0	0.0809231867781842
            1182	0.0606524478467524	0	0	0.0606524478467524
            1183	0.142175934884839	0	0	0.142175934884839
            1184	0.190535015883832	0	0	0.190535015883832
            1185	0.376721142545183	0	0	0.376721142545183
            1186	0.281039126181166	0	0	0.281039126181166
            1187	0.457907935484970	0	0	0.457907935484970
            1188	0.613781606222243	0	0	0.613781606222243
            1189	0.822940089594963	0	0	0.822940089594963
            1190	0.613971823059841	0	0	0.613971823059841
            1191	0.677995264599961	0	0	0.677995264599961
            1192	0.908730420862152	0	0	0.908730420862152
            1193	0.830165946767061	0	0	0.830165946767061
            1194	0.619393154246033	0	0	0.619393154246033
            1195	0.478006422071870	0	0	0.478006422071870
            1196	0.640652410691288	0	0	0.640652410691288
            1197	0.455550008781949	0	0	0.455550008781949
            1198	0.339901988229740	0	0	0.339901988229740
            1199	0.283691548385499	0	0	0.283691548385499
            1200	0.380222524874520	0	0	0.380222524874520
            1201	0.455546210185846	0	0	0.455546210185846
            1202	0.339875140256020	0	0	0.339875140256020
            1203	0.477959208312821	0	0	0.477959208312821
            1204	0.640645133442157	0	0	0.640645133442157
            1205	0.830155490407122	0	0	0.830155490407122
            1206	0.619337042288007	0	0	0.619337042288007
            1207	0.677942729443061	0	0	0.677942729443061
            1208	0.908716686440818	0	0	0.908716686440818
            1209	0.822922554696522	0	0	0.822922554696522
            1210	0.613930920741938	0	0	0.613930920741938
            1211	0.457876028938820	0	0	0.457876028938820
            1212	0.613756483578898	0	0	0.613756483578898
            1213	0.376652051329979	0	0	0.376652051329979
            1214	0.280978343050236	0	0	0.280978343050236
            1215	0.141876230643514	0	0	0.141876230643514
            1216	0.190185865525938	0	0	0.190185865525938
            1217	0.0791502187346933	0	0	0.0791502187346933
            1218	0.0591488616399993	0	0	0.0591488616399993
            1219	0.0500000000000000	0	0	0.0500000000000000
            1220	0.0500000000000000	0	0	0.0500000000000000
            1221	0.0500000000000000	0	0	0.0500000000000000
            1222	0.0500000000000000	0	0	0.0500000000000000
            1223	0.0500000000000000	0	0	0.0500000000000000
            1224	0.0500000000000000	0	0	0.0500000000000000
            1225	0.0500000000000000	0	0	0.0500000000000000
            1226	0.0591718082103400	0	0	0.0591718082103400
            1227	0.142060880870281	0	0	0.142060880870281
            1228	0.0873675697234437	0	0	0.0873675697234437
            1229	0.173888632955907	0	0	0.173888632955907
            1230	0.282018493538424	0	0	0.282018493538424
            1231	0.462608724230332	0	0	0.462608724230332
            1232	0.287417889935360	0	0	0.287417889935360
            1233	0.400187843340956	0	0	0.400187843340956
            1234	0.631677101354358	0	0	0.631677101354358
            1235	0.732891713212133	0	0	0.732891713212133
            1236	0.486474671589488	0	0	0.486474671589488
            1237	0.554575214070496	0	0	0.554575214070496
            1238	0.759612454478278	0	0	0.759612454478278
            1239	0.772086060655468	0	0	0.772086060655468
            1240	0.648965042065232	0	0	0.648965042065232
            1241	0.802676782061594	0	0	0.802676782061594
            1242	0.842891294328642	0	0	0.842891294328642
            1243	0.950000000000000	0	0	0.950000000000000
            1244	0.950000000000000	0	0	0.950000000000000
            1245	0.950000000000000	0	0	0.950000000000000
            1246	0.950000000000000	0	0	0.950000000000000
            1247	0.950000000000000	0	0	0.950000000000000
            1248	0.950000000000000	0	0	0.950000000000000
            1249	0.837917246623500	0	0	0.837917246623500
            1250	0.900372338206990	0	0	0.900372338206990
            1251	0.631713610412352	0	0	0.631713610412352
            1252	0.562752015426341	0	0	0.562752015426341
            1253	0.314443012855653	0	0	0.314443012855653
            1254	0.368498312589345	0	0	0.368498312589345
            1255	0.178343382163130	0	0	0.178343382163130
            1256	0.146007369607714	0	0	0.146007369607714
            1257	0.0563178144347014	0	0	0.0563178144347014
            1258	0.0716097663216013	0	0	0.0716097663216013
            1259	0.0500000000000000	0	0	0.0500000000000000
            1260	0.0500000000000000	0	0	0.0500000000000000
            1261	0.0500000000000000	0	0	0.0500000000000000
            1262	0.0500000000000000	0	0	0.0500000000000000
            1263	0.0500000000000000	0	0	0.0500000000000000
            1264	0.0500000000000000	0	0	0.0500000000000000
            1265	0.0500000000000000	0	0	0.0500000000000000
            1266	0.0500000000000000	0	0	0.0500000000000000
            1267	0.0500000000000000	0	0	0.0500000000000000
            1268	0.0500000000000000	0	0	0.0500000000000000
            1269	0.0500000000000000	0	0	0.0500000000000000
            1270	0.0500000000000000	0	0	0.0500000000000000
            1271	0.0500000000000000	0	0	0.0500000000000000
            1272	0.0500000000000000	0	0	0.0500000000000000
            1273	0.0500000000000000	0	0	0.0500000000000000
            1274	0.0500000000000000	0	0	0.0500000000000000
            1275	0.0500000000000000	0	0	0.0500000000000000
            1276	0.0500000000000000	0	0	0.0500000000000000
            1277	0.0500000000000000	0	0	0.0500000000000000
            1278	0.0500000000000000	0	0	0.0500000000000000
            1279	0.0500000000000000	0	0	0.0500000000000000
            1280	0.0500000000000000	0	0	0.0500000000000000
            1281	0.842693884977158	0	0	0.842693884977158
            1282	0.919310278627720	0	0	0.919310278627720
            1283	0.950000000000000	0	0	0.950000000000000
            1284	0.950000000000000	0	0	0.950000000000000
            1285	0.950000000000000	0	0	0.950000000000000
            1286	0.950000000000000	0	0	0.950000000000000
            1287	0.849062165649411	0	0	0.849062165649411
            1288	0.894644252135583	0	0	0.894644252135583
            1289	0.834437341160397	0	0	0.834437341160397
            1290	0.696205116435428	0	0	0.696205116435428
            1291	0.586217904409247	0	0	0.586217904409247
            1292	0.803789850836583	0	0	0.803789850836583
            1293	0.755574410135824	0	0	0.755574410135824
            1294	0.502276758191438	0	0	0.502276758191438
            1295	0.406333473865232	0	0	0.406333473865232
            1296	0.640643873556393	0	0	0.640643873556393
            1297	0.465412134727122	0	0	0.465412134727122
            1298	0.289317880296415	0	0	0.289317880296415
            1299	0.174364186243505	0	0	0.174364186243505
            1300	0.282738001627410	0	0	0.282738001627410
            1301	0.142315492498069	0	0	0.142315492498069
            1302	0.0874824234749989	0	0	0.0874824234749989
            1303	0.0500000000000000	0	0	0.0500000000000000
            1304	0.0597470512403920	0	0	0.0597470512403920
            1305	0.0500000000000000	0	0	0.0500000000000000
            1306	0.0500000000000000	0	0	0.0500000000000000
            1307	0.0500000000000000	0	0	0.0500000000000000
            1308	0.0500000000000000	0	0	0.0500000000000000
            1309	0.0500000000000000	0	0	0.0500000000000000
            1310	0.0500000000000000	0	0	0.0500000000000000
            1311	0.0500000000000000	0	0	0.0500000000000000
            1312	0.0873198542577423	0	0	0.0873198542577423
            1313	0.172541637144458	0	0	0.172541637144458
            1314	0.0874715230283261	0	0	0.0874715230283261
            1315	0.142657809898740	0	0	0.142657809898740
            1316	0.281146744193454	0	0	0.281146744193454
            1317	0.377010012471161	0	0	0.377010012471161
            1318	0.191524708094005	0	0	0.191524708094005
            1319	0.211843848508462	0	0	0.211843848508462
            1320	0.416383515959216	0	0	0.416383515959216
            1321	0.380462747935907	0	0	0.380462747935907
            1322	0.193942972351784	0	0	0.193942972351784
            1323	0.150060912284419	0	0	0.150060912284419
            1324	0.293679666076147	0	0	0.293679666076147
            1325	0.208863540258662	0	0	0.208863540258662
            1326	0.106921663898095	0	0	0.106921663898095
            1327	0.0890424602984019	0	0	0.0890424602984019
            1328	0.174286518762751	0	0	0.174286518762751
            1329	0.208692139522053	0	0	0.208692139522053
            1330	0.106020098000341	0	0	0.106020098000341
            1331	0.148478088284478	0	0	0.148478088284478
            1332	0.293378740029589	0	0	0.293378740029589
            1333	0.380106613421932	0	0	0.380106613421932
            1334	0.192069903047045	0	0	0.192069903047045
            1335	0.210109094638023	0	0	0.210109094638023
            1336	0.416053591985692	0	0	0.416053591985692
            1337	0.376761770560471	0	0	0.376761770560471
            1338	0.190221951740461	0	0	0.190221951740461
            1339	0.141854617263725	0	0	0.141854617263725
            1340	0.280989996975846	0	0	0.280989996975846
            1341	0.172432413931411	0	0	0.172432413931411
            1342	0.0870505198021727	0	0	0.0870505198021727
            1343	0.0500000000000000	0	0	0.0500000000000000
            1344	0.0870915653432057	0	0	0.0870915653432057
            1345	0.0500000000000000	0	0	0.0500000000000000
            1346	0.0500000000000000	0	0	0.0500000000000000
            1347	0.0500000000000000	0	0	0.0500000000000000
            1348	0.0500000000000000	0	0	0.0500000000000000
            1349	0.0500000000000000	0	0	0.0500000000000000
            1350	0.0500000000000000	0	0	0.0500000000000000
            1351	0.0500000000000000	0	0	0.0500000000000000
            1352	0.0500000000000000	0	0	0.0500000000000000
            1353	0.0791117132508334	0	0	0.0791117132508334
            1354	0.0871281657597233	0	0	0.0871281657597233
            1355	0.209714462774774	0	0	0.209714462774774
            1356	0.190262369559574	0	0	0.190262369559574
            1357	0.377344679168920	0	0	0.377344679168920
            1358	0.415778394209706	0	0	0.415778394209706
            1359	0.679303510646882	0	0	0.679303510646882
            1360	0.617183821414025	0	0	0.617183821414025
            1361	0.836387918838445	0	0	0.836387918838445
            1362	0.917924870590799	0	0	0.917924870590799
            1363	0.950000000000000	0	0	0.950000000000000
            1364	0.950000000000000	0	0	0.950000000000000
            1365	0.945461674270938	0	0	0.945461674270938
            1366	0.950000000000000	0	0	0.950000000000000
            1367	0.922447225429840	0	0	0.922447225429840
            1368	0.890835982988891	0	0	0.890835982988891
            1369	0.895047482750403	0	0	0.895047482750403
            1370	0.889434938848769	0	0	0.889434938848769
            1371	0.950000000000000	0	0	0.950000000000000
            1372	0.950000000000000	0	0	0.950000000000000
            1373	0.950000000000000	0	0	0.950000000000000
            1374	0.950000000000000	0	0	0.950000000000000
            1375	0.950000000000000	0	0	0.950000000000000
            1376	0.950000000000000	0	0	0.950000000000000
            1377	0.950000000000000	0	0	0.950000000000000
            1378	0.950000000000000	0	0	0.950000000000000
            1379	0.718175940541122	0	0	0.718175940541122
            1380	0.703271279913444	0	0	0.703271279913444
            1381	0.422179508956860	0	0	0.422179508956860
            1382	0.439500340703347	0	0	0.439500340703347
            1383	0.224820257975447	0	0	0.224820257975447
            1384	0.210156419160674	0	0	0.210156419160674
            1385	0.0871091130689944	0	0	0.0871091130689944
            1386	0.0979862338146197	0	0	0.0979862338146197
            1387	0.0500000000000000	0	0	0.0500000000000000
            1388	0.0500000000000000	0	0	0.0500000000000000
            1389	0.0500000000000000	0	0	0.0500000000000000
            1390	0.0500000000000000	0	0	0.0500000000000000
            1391	0.0500000000000000	0	0	0.0500000000000000
            1392	0.0500000000000000	0	0	0.0500000000000000
            1393	0.0500000000000000	0	0	0.0500000000000000
            1394	0.0500000000000000	0	0	0.0500000000000000
            1395	0.0500000000000000	0	0	0.0500000000000000
            1396	0.0500000000000000	0	0	0.0500000000000000
            1397	0.0500000000000000	0	0	0.0500000000000000
            1398	0.0500000000000000	0	0	0.0500000000000000
            1399	0.0500000000000000	0	0	0.0500000000000000
            1400	0.0500000000000000	0	0	0.0500000000000000
            1401	0.0500000000000000	0	0	0.0500000000000000
            1402	0.0500000000000000	0	0	0.0500000000000000
            1403	0.0500000000000000	0	0	0.0500000000000000
            1404	0.0500000000000000	0	0	0.0500000000000000
            1405	0.0500000000000000	0	0	0.0500000000000000
            1406	0.0500000000000000	0	0	0.0500000000000000
            1407	0.0500000000000000	0	0	0.0500000000000000
            1408	0.0500000000000000	0	0	0.0500000000000000
            1409	0.828847977304331	0	0	0.828847977304331
            1410	0.616516156440097	0	0	0.616516156440097
            1411	0.684542375775862	0	0	0.684542375775862
            1412	0.925256902153044	0	0	0.925256902153044
            1413	0.866571886951586	0	0	0.866571886951586
            1414	0.632291481980188	0	0	0.632291481980188
            1415	0.494307751594017	0	0	0.494307751594017
            1416	0.700770670276232	0	0	0.700770670276232
            1417	0.518472442866586	0	0	0.518472442866586
            1418	0.340104726636059	0	0	0.340104726636059
            1419	0.220922936144516	0	0	0.220922936144516
            1420	0.380242378684269	0	0	0.380242378684269
            1421	0.288563796691894	0	0	0.288563796691894
            1422	0.146008431178995	0	0	0.146008431178995
            1423	0.0990881833549943	0	0	0.0990881833549943
            1424	0.217077569016427	0	0	0.217077569016427
            1425	0.149408197246041	0	0	0.149408197246041
            1426	0.0645561325400273	0	0	0.0645561325400273
            1427	0.0500000000000000	0	0	0.0500000000000000
            1428	0.0887947578450023	0	0	0.0887947578450023
            1429	0.0500000000000000	0	0	0.0500000000000000
            1430	0.0500000000000000	0	0	0.0500000000000000
            1431	0.0500000000000000	0	0	0.0500000000000000
            1432	0.0500000000000000	0	0	0.0500000000000000
            1433	0.0500000000000000	0	0	0.0500000000000000
            1434	0.0500000000000000	0	0	0.0500000000000000
            1435	0.0500000000000000	0	0	0.0500000000000000
            1436	0.0500000000000000	0	0	0.0500000000000000
            1437	0.0500000000000000	0	0	0.0500000000000000
            1438	0.0500000000000000	0	0	0.0500000000000000
            1439	0.0500000000000000	0	0	0.0500000000000000
            1440	0.0500000000000000	0	0	0.0500000000000000
            1441	0.0500000000000000	0	0	0.0500000000000000
            1442	0.0500000000000000	0	0	0.0500000000000000
            1443	0.0500000000000000	0	0	0.0500000000000000
            1444	0.0624254872250737	0	0	0.0624254872250737
            1445	0.0847923322637202	0	0	0.0847923322637202
            1446	0.0500000000000000	0	0	0.0500000000000000
            1447	0.0585008533549907	0	0	0.0585008533549907
            1448	0.0951867435153072	0	0	0.0951867435153072
            1449	0.0888218089728071	0	0	0.0888218089728071
            1450	0.0604092344774486	0	0	0.0604092344774486
            1451	0.0536146719046674	0	0	0.0536146719046674
            1452	0.0703625782237588	0	0	0.0703625782237588
            1453	0.0511198594133398	0	0	0.0511198594133398
            1454	0.0500000000000000	0	0	0.0500000000000000
            1455	0.0500000000000000	0	0	0.0500000000000000
            1456	0.0500000000000000	0	0	0.0500000000000000
            1457	0.0500000000000000	0	0	0.0500000000000000
            1458	0.0500000000000000	0	0	0.0500000000000000
            1459	0.0500000000000000	0	0	0.0500000000000000
            1460	0.0635139142902376	0	0	0.0635139142902376
            1461	0.0807173062062752	0	0	0.0807173062062752
            1462	0.0500000000000000	0	0	0.0500000000000000
            1463	0.0500000000000000	0	0	0.0500000000000000
            1464	0.0876808108181858	0	0	0.0876808108181858
            1465	0.0791565255440016	0	0	0.0791565255440016
            1466	0.0500000000000000	0	0	0.0500000000000000
            1467	0.0500000000000000	0	0	0.0500000000000000
            1468	0.0589592675245787	0	0	0.0589592675245787
            1469	0.0500000000000000	0	0	0.0500000000000000
            1470	0.0500000000000000	0	0	0.0500000000000000
            1471	0.0500000000000000	0	0	0.0500000000000000
            1472	0.0500000000000000	0	0	0.0500000000000000
            1473	0.0500000000000000	0	0	0.0500000000000000
            1474	0.0500000000000000	0	0	0.0500000000000000
            1475	0.0500000000000000	0	0	0.0500000000000000
            1476	0.0500000000000000	0	0	0.0500000000000000
            1477	0.0500000000000000	0	0	0.0500000000000000
            1478	0.0500000000000000	0	0	0.0500000000000000
            1479	0.0500000000000000	0	0	0.0500000000000000
            1480	0.0500000000000000	0	0	0.0500000000000000
            1481	0.0789837829460725	0	0	0.0789837829460725
            1482	0.0589256162196319	0	0	0.0589256162196319
            1483	0.141931154619995	0	0	0.141931154619995
            1484	0.190188141229885	0	0	0.190188141229885
            1485	0.377026279337230	0	0	0.377026279337230
            1486	0.281395983688340	0	0	0.281395983688340
            1487	0.459614016463188	0	0	0.459614016463188
            1488	0.615756416685745	0	0	0.615756416685745
            1489	0.831189486073429	0	0	0.831189486073429
            1490	0.620382039451861	0	0	0.620382039451861
            1491	0.698783509244751	0	0	0.698783509244751
            1492	0.936759971468373	0	0	0.936759971468373
            1493	0.908084074644967	0	0	0.908084074644967
            1494	0.676080043664864	0	0	0.676080043664864
            1495	0.605094315718141	0	0	0.605094315718141
            1496	0.816512097714701	0	0	0.816512097714701
            1497	0.773712938378857	0	0	0.773712938378857
            1498	0.569497254096400	0	0	0.569497254096400
            1499	0.609272424603672	0	0	0.609272424603672
            1500	0.831249398370034	0	0	0.831249398370034
            1501	0.932283367573391	0	0	0.932283367573391
            1502	0.684315353332035	0	0	0.684315353332035
            1503	0.712305640531978	0	0	0.712305640531978
            1504	0.950000000000000	0	0	0.950000000000000
            1505	0.856532559418870	0	0	0.856532559418870
            1506	0.643396109903486	0	0	0.643396109903486
            1507	0.498806917960264	0	0	0.498806917960264
            1508	0.640131341971987	0	0	0.640131341971987
            1509	0.403286765585318	0	0	0.403286765585318
            1510	0.341871825975034	0	0	0.341871825975034
            1511	0.221540814227886	0	0	0.221540814227886
            1512	0.219712954674719	0	0	0.219712954674719
            1513	0.109560724049023	0	0	0.109560724049023
            1514	0.146188565017626	0	0	0.146188565017626
            1515	0.0991374851952626	0	0	0.0991374851952626
            1516	0.0541772153952575	0	0	0.0541772153952575
            1517	0.0500000000000000	0	0	0.0500000000000000
            1518	0.0645930330196613	0	0	0.0645930330196613
            1519	0.0500000000000000	0	0	0.0500000000000000
            1520	0.0500000000000000	0	0	0.0500000000000000
            1521	0.0500000000000000	0	0	0.0500000000000000
            1522	0.0500000000000000	0	0	0.0500000000000000
            1523	0.0500000000000000	0	0	0.0500000000000000
            1524	0.0500000000000000	0	0	0.0500000000000000
            1525	0.0500000000000000	0	0	0.0500000000000000
            1526	0.0500000000000000	0	0	0.0500000000000000
            1527	0.0500000000000000	0	0	0.0500000000000000
            1528	0.0500000000000000	0	0	0.0500000000000000
            1529	0.0500000000000000	0	0	0.0500000000000000
            1530	0.0500000000000000	0	0	0.0500000000000000
            1531	0.0500000000000000	0	0	0.0500000000000000
            1532	0.0500000000000000	0	0	0.0500000000000000
            1533	0.0500000000000000	0	0	0.0500000000000000
            1534	0.0500000000000000	0	0	0.0500000000000000
            1535	0.0500000000000000	0	0	0.0500000000000000
            1536	0.0500000000000000	0	0	0.0500000000000000
            1537	0.378793487913849	0	0	0.378793487913849
            1538	0.196009398388013	0	0	0.196009398388013
            1539	0.217914094543965	0	0	0.217914094543965
            1540	0.419575258534704	0	0	0.419575258534704
            1541	0.384765493665121	0	0	0.384765493665121
            1542	0.200212999325082	0	0	0.200212999325082
            1543	0.152874454347183	0	0	0.152874454347183
            1544	0.295058827140732	0	0	0.295058827140732
            1545	0.193843696862199	0	0	0.193843696862199
            1546	0.0983426127889171	0	0	0.0983426127889171
            1547	0.0549825739980459	0	0	0.0549825739980459
            1548	0.114803357507467	0	0	0.114803357507467
            1549	0.0662757319023152	0	0	0.0662757319023152
            1550	0.0500000000000000	0	0	0.0500000000000000
            1551	0.0500000000000000	0	0	0.0500000000000000
            1552	0.0500000000000000	0	0	0.0500000000000000
            1553	0.0500000000000000	0	0	0.0500000000000000
            1554	0.0500000000000000	0	0	0.0500000000000000
            1555	0.0500000000000000	0	0	0.0500000000000000
            1556	0.0500000000000000	0	0	0.0500000000000000
            1557	0.0500000000000000	0	0	0.0500000000000000
            1558	0.0500000000000000	0	0	0.0500000000000000
            1559	0.0500000000000000	0	0	0.0500000000000000
            1560	0.0500000000000000	0	0	0.0500000000000000
            1561	0.0500000000000000	0	0	0.0500000000000000
            1562	0.0500000000000000	0	0	0.0500000000000000
            1563	0.0500000000000000	0	0	0.0500000000000000
            1564	0.0500000000000000	0	0	0.0500000000000000
            1565	0.0500000000000000	0	0	0.0500000000000000
            1566	0.0500000000000000	0	0	0.0500000000000000
            1567	0.0500000000000000	0	0	0.0500000000000000
            1568	0.0500000000000000	0	0	0.0500000000000000
            1569	0.0500000000000000	0	0	0.0500000000000000
            1570	0.0500000000000000	0	0	0.0500000000000000
            1571	0.0909727275022944	0	0	0.0909727275022944
            1572	0.0500000000000000	0	0	0.0500000000000000
            1573	0.0690625878093345	0	0	0.0690625878093345
            1574	0.151181873479933	0	0	0.151181873479933
            1575	0.210299937225882	0	0	0.210299937225882
            1576	0.0935424699820199	0	0	0.0935424699820199
            1577	0.107029129441457	0	0	0.107029129441457
            1578	0.247524950488628	0	0	0.247524950488628
            1579	0.250367700619668	0	0	0.250367700619668
            1580	0.104956509197479	0	0	0.104956509197479
            1581	0.0900518094465889	0	0	0.0900518094465889
            1582	0.221740632113056	0	0	0.221740632113056
            1583	0.174886136865651	0	0	0.174886136865651
            1584	0.0693724581942503	0	0	0.0693724581942503
            1585	0.0500000000000000	0	0	0.0500000000000000
            1586	0.123846945680540	0	0	0.123846945680540
            1587	0.0785009432501049	0	0	0.0785009432501049
            1588	0.0500000000000000	0	0	0.0500000000000000
            1589	0.0500000000000000	0	0	0.0500000000000000
            1590	0.0500000000000000	0	0	0.0500000000000000
            1591	0.0500000000000000	0	0	0.0500000000000000
            1592	0.0500000000000000	0	0	0.0500000000000000
            1593	0.0500000000000000	0	0	0.0500000000000000
            1594	0.0500000000000000	0	0	0.0500000000000000
            1595	0.0500000000000000	0	0	0.0500000000000000
            1596	0.0500000000000000	0	0	0.0500000000000000
            1597	0.0500000000000000	0	0	0.0500000000000000
            1598	0.0500000000000000	0	0	0.0500000000000000
            1599	0.0500000000000000	0	0	0.0500000000000000
            1600	0.0500000000000000	0	0	0.0500000000000000
            1601	0.0500000000000000	0	0	0.0500000000000000
            1602	0.0500000000000000	0	0	0.0500000000000000
            1603	0.0500000000000000	0	0	0.0500000000000000
            1604	0.0500000000000000	0	0	0.0500000000000000
            1605	0.0500000000000000	0	0	0.0500000000000000
            1606	0.0500000000000000	0	0	0.0500000000000000
            1607	0.0500000000000000	0	0	0.0500000000000000
            1608	0.0500000000000000	0	0	0.0500000000000000
            1609	0.0500000000000000	0	0	0.0500000000000000
            1610	0.0500000000000000	0	0	0.0500000000000000
            1611	0.0500000000000000	0	0	0.0500000000000000
            1612	0.0873117141376390	0	0	0.0873117141376390
            1613	0.173304824112184	0	0	0.173304824112184
            1614	0.0897036633546457	0	0	0.0897036633546457
            1615	0.148469527748699	0	0	0.148469527748699
            1616	0.283531882792293	0	0	0.283531882792293
            1617	0.383597803669905	0	0	0.383597803669905
            1618	0.204216047878701	0	0	0.204216047878701
            1619	0.235993786024467	0	0	0.235993786024467
            1620	0.433402866013892	0	0	0.433402866013892
            1621	0.420737339293009	0	0	0.420737339293009
            1622	0.235417498176370	0	0	0.235417498176370
            1623	0.215801092651903	0	0	0.215801092651903
            1624	0.377233001231021	0	0	0.377233001231021
            1625	0.354129781050243	0	0	0.354129781050243
            1626	0.202089882748015	0	0	0.202089882748015
            1627	0.209684136818074	0	0	0.209684136818074
            1628	0.377155151607051	0	0	0.377155151607051
            1629	0.424375420896349	0	0	0.424375420896349
            1630	0.233409835074148	0	0	0.233409835074148
            1631	0.258548946045932	0	0	0.258548946045932
            1632	0.449327420452083	0	0	0.449327420452083
            1633	0.426042022752128	0	0	0.426042022752128
            1634	0.281734453617472	0	0	0.281734453617472
            1635	0.315693124598501	0	0	0.315693124598501
            1636	0.369853284960081	0	0	0.369853284960081
            1637	0.314879296964439	0	0	0.314879296964439
            1638	0.368406901921760	0	0	0.368406901921760
            1639	0.420657872577646	0	0	0.420657872577646
            1640	0.277359346635028	0	0	0.277359346635028
            1641	0.245832569453047	0	0	0.245832569453047
            1642	0.433574199760684	0	0	0.433574199760684
            1643	0.383088291129672	0	0	0.383088291129672
            1644	0.202528258096355	0	0	0.202528258096355
            1645	0.145470397085101	0	0	0.145470397085101
            1646	0.283442096221327	0	0	0.283442096221327
            1647	0.174577273870292	0	0	0.174577273870292
            1648	0.0882453254965428	0	0	0.0882453254965428
            1649	0.0500000000000000	0	0	0.0500000000000000
            1650	0.0905837154716453	0	0	0.0905837154716453
            1651	0.0500000000000000	0	0	0.0500000000000000
            1652	0.0500000000000000	0	0	0.0500000000000000
            1653	0.0500000000000000	0	0	0.0500000000000000
            1654	0.0500000000000000	0	0	0.0500000000000000
            1655	0.0500000000000000	0	0	0.0500000000000000
            1656	0.0500000000000000	0	0	0.0500000000000000
            1657	0.0500000000000000	0	0	0.0500000000000000
            1658	0.0500000000000000	0	0	0.0500000000000000
            1659	0.0500000000000000	0	0	0.0500000000000000
            1660	0.0500000000000000	0	0	0.0500000000000000
            1661	0.0500000000000000	0	0	0.0500000000000000
            1662	0.0500000000000000	0	0	0.0500000000000000
            1663	0.0500000000000000	0	0	0.0500000000000000
            1664	0.0500000000000000	0	0	0.0500000000000000
            1665	0.0991135776925834	0	0	0.0991135776925834
            1666	0.0858673791872043	0	0	0.0858673791872043
            1667	0.108702239704177	0	0	0.108702239704177
            1668	0.114169730534565	0	0	0.114169730534565
            1669	0.109074713090502	0	0	0.109074713090502
            1670	0.114081591221354	0	0	0.114081591221354
            1671	0.0992194896152273	0	0	0.0992194896152273
            1672	0.0866596393549081	0	0	0.0866596393549081
            1673	0.0575685593803108	0	0	0.0575685593803108
            1674	0.0715181855402558	0	0	0.0715181855402558
            1675	0.0500000000000000	0	0	0.0500000000000000
            1676	0.0500000000000000	0	0	0.0500000000000000
            1677	0.0500000000000000	0	0	0.0500000000000000
            1678	0.0500000000000000	0	0	0.0500000000000000
            1679	0.0500000000000000	0	0	0.0500000000000000
            1680	0.0500000000000000	0	0	0.0500000000000000
            1681	0.0500000000000000	0	0	0.0500000000000000
            1682	0.0500000000000000	0	0	0.0500000000000000
            1683	0.0500000000000000	0	0	0.0500000000000000
            1684	0.0500000000000000	0	0	0.0500000000000000
            1685	0.0500000000000000	0	0	0.0500000000000000
            1686	0.0500000000000000	0	0	0.0500000000000000
            1687	0.0500000000000000	0	0	0.0500000000000000
            1688	0.0500000000000000	0	0	0.0500000000000000
            1689	0.0500000000000000	0	0	0.0500000000000000
            1690	0.0500000000000000	0	0	0.0500000000000000
            1691	0.0500000000000000	0	0	0.0500000000000000
            1692	0.0500000000000000	0	0	0.0500000000000000
            1693	0.0500000000000000	0	0	0.0500000000000000
            1694	0.0500000000000000	0	0	0.0500000000000000
            1695	0.0618965591275142	0	0	0.0618965591275142
            1696	0.0500000000000000	0	0	0.0500000000000000
            1697	0.0892988801919474	0	0	0.0892988801919474
            1698	0.147228711792111	0	0	0.147228711792111
            1699	0.294730642513222	0	0	0.294730642513222
            1700	0.179162027575041	0	0	0.179162027575041
            1701	0.299697920709589	0	0	0.299697920709589
            1702	0.494860643336777	0	0	0.494860643336777
            1703	0.701191697330262	0	0	0.701191697330262
            1704	0.420986726448628	0	0	0.420986726448628
            1705	0.502643536485540	0	0	0.502643536485540
            1706	0.849279092456589	0	0	0.849279092456589
            1707	0.894742182110068	0	0	0.894742182110068
            1708	0.518703864464075	0	0	0.518703864464075
            1709	0.471259985801386	0	0	0.471259985801386
            1710	0.834481897956948	0	0	0.834481897956948
            1711	0.696227928466483	0	0	0.696227928466483
            1712	0.382187576410916	0	0	0.382187576410916
            1713	0.277332109589576	0	0	0.277332109589576
            1714	0.518463147144129	0	0	0.518463147144129
            1715	0.339897048240696	0	0	0.339897048240696
            1716	0.178230949269748	0	0	0.178230949269748
            1717	0.0997189534899161	0	0	0.0997189534899161
            1718	0.192589963931774	0	0	0.192589963931774
            1719	0.0927476731577740	0	0	0.0927476731577740
            1720	0.0500000000000000	0	0	0.0500000000000000
            1721	0.0500000000000000	0	0	0.0500000000000000
            1722	0.0500000000000000	0	0	0.0500000000000000
            1723	0.0500000000000000	0	0	0.0500000000000000
            1724	0.0500000000000000	0	0	0.0500000000000000
            1725	0.0500000000000000	0	0	0.0500000000000000
            1726	0.0500000000000000	0	0	0.0500000000000000
            1727	0.0500000000000000	0	0	0.0500000000000000
            1728	0.0500000000000000	0	0	0.0500000000000000
            1729	0.0500000000000000	0	0	0.0500000000000000
            1730	0.0500000000000000	0	0	0.0500000000000000
            1731	0.0500000000000000	0	0	0.0500000000000000
            1732	0.0500000000000000	0	0	0.0500000000000000
            1733	0.0500000000000000	0	0	0.0500000000000000
            1734	0.0500000000000000	0	0	0.0500000000000000
            1735	0.0500000000000000	0	0	0.0500000000000000
            1736	0.0500000000000000	0	0	0.0500000000000000
            1737	0.0500000000000000	0	0	0.0500000000000000
            1738	0.0500000000000000	0	0	0.0500000000000000
            1739	0.0500000000000000	0	0	0.0500000000000000
            1740	0.0500000000000000	0	0	0.0500000000000000
            1741	0.0500000000000000	0	0	0.0500000000000000
            1742	0.0500000000000000	0	0	0.0500000000000000
            1743	0.0641152371287546	0	0	0.0641152371287546
            1744	0.0773064162486711	0	0	0.0773064162486711
            1745	0.115811482281852	0	0	0.115811482281852
            1746	0.114253912278916	0	0	0.114253912278916
            1747	0.172452087361944	0	0	0.172452087361944
            1748	0.148529209749933	0	0	0.148529209749933
            1749	0.165507351253387	0	0	0.165507351253387
            1750	0.219714325409429	0	0	0.219714325409429
            1751	0.236294900435781	0	0	0.236294900435781
            1752	0.164447882888935	0	0	0.164447882888935
            1753	0.152568492985205	0	0	0.152568492985205
            1754	0.216554507043283	0	0	0.216554507043283
            1755	0.175759923101862	0	0	0.175759923101862
            1756	0.141734296569790	0	0	0.141734296569790
            1757	0.142830772135861	0	0	0.142830772135861
            1758	0.143559798116043	0	0	0.143559798116043
            1759	0.153497112023200	0	0	0.153497112023200
            1760	0.166815363116070	0	0	0.166815363116070
            1761	0.229112276027423	0	0	0.229112276027423
            1762	0.235570695218939	0	0	0.235570695218939
            1763	0.403016278589148	0	0	0.403016278589148
            1764	0.342637540969396	0	0	0.342637540969396
            1765	0.494756796382811	0	0	0.494756796382811
            1766	0.627820163430623	0	0	0.627820163430623
            1767	0.829384162172272	0	0	0.829384162172272
            1768	0.632159889060537	0	0	0.632159889060537
            1769	0.684492137935503	0	0	0.684492137935503
            1770	0.910352581791707	0	0	0.910352581791707
            1771	0.826061085851289	0	0	0.826061085851289
            1772	0.617157202546039	0	0	0.617157202546039
            1773	0.461148719192461	0	0	0.461148719192461
            1774	0.621505815134412	0	0	0.621505815134412
            1775	0.394905241232778	0	0	0.394905241232778
            1776	0.287362987801580	0	0	0.287362987801580
            1777	0.154206776701655	0	0	0.154206776701655
            1778	0.226261734309293	0	0	0.226261734309293
            1779	0.137835184634926	0	0	0.137835184634926
            1780	0.0790138955647173	0	0	0.0790138955647173
            1781	0.0500000000000000	0	0	0.0500000000000000
            1782	0.105903819971215	0	0	0.105903819971215
            1783	0.0946140686124669	0	0	0.0946140686124669
            1784	0.0500000000000000	0	0	0.0500000000000000
            1785	0.0500000000000000	0	0	0.0500000000000000
            1786	0.0806912198662532	0	0	0.0806912198662532
            1787	0.0592259578460571	0	0	0.0592259578460571
            1788	0.0500000000000000	0	0	0.0500000000000000
            1789	0.0500000000000000	0	0	0.0500000000000000
            1790	0.0500000000000000	0	0	0.0500000000000000
            1791	0.0500000000000000	0	0	0.0500000000000000
            1792	0.0500000000000000	0	0	0.0500000000000000
            1793	0.149408251323025	0	0	0.149408251323025
            1794	0.282717423846731	0	0	0.282717423846731
            1795	0.378550678805720	0	0	0.378550678805720
            1796	0.198470815136843	0	0	0.198470815136843
            1797	0.217215233933025	0	0	0.217215233933025
            1798	0.417059055265666	0	0	0.417059055265666
            1799	0.378161552129225	0	0	0.378161552129225
            1800	0.195842411943719	0	0	0.195842411943719
            1801	0.145506183303435	0	0	0.145506183303435
            1802	0.282500183961690	0	0	0.282500183961690
            1803	0.174601641729141	0	0	0.174601641729141
            1804	0.0892501162674293	0	0	0.0892501162674293
            1805	0.0500000000000000	0	0	0.0500000000000000
            1806	0.0907628565319331	0	0	0.0907628565319331
            1807	0.0500000000000000	0	0	0.0500000000000000
            1808	0.0500000000000000	0	0	0.0500000000000000
            1809	0.0500000000000000	0	0	0.0500000000000000
            1810	0.0500000000000000	0	0	0.0500000000000000
            1811	0.0500000000000000	0	0	0.0500000000000000
            1812	0.0500000000000000	0	0	0.0500000000000000
            1813	0.0500000000000000	0	0	0.0500000000000000
            1814	0.0500000000000000	0	0	0.0500000000000000
            1815	0.0500000000000000	0	0	0.0500000000000000
            1816	0.0500000000000000	0	0	0.0500000000000000
            1817	0.0500000000000000	0	0	0.0500000000000000
            1818	0.0500000000000000	0	0	0.0500000000000000
            1819	0.0500000000000000	0	0	0.0500000000000000
            1820	0.0500000000000000	0	0	0.0500000000000000
            1821	0.0500000000000000	0	0	0.0500000000000000
            1822	0.0653283967538326	0	0	0.0653283967538326
            1823	0.124418747076591	0	0	0.124418747076591
            1824	0.0899645851337002	0	0	0.0899645851337002
            1825	0.204157540589253	0	0	0.204157540589253
            1826	0.247688787320029	0	0	0.247688787320029
            1827	0.464083350550270	0	0	0.464083350550270
            1828	0.402784412396223	0	0	0.402784412396223
            1829	0.676217184928639	0	0	0.676217184928639
            1830	0.768649673322300	0	0	0.768649673322300
            1831	0.950000000000000	0	0	0.950000000000000
            1832	0.950000000000000	0	0	0.950000000000000
            1833	0.950000000000000	0	0	0.950000000000000
            1834	0.950000000000000	0	0	0.950000000000000
            1835	0.950000000000000	0	0	0.950000000000000
            1836	0.950000000000000	0	0	0.950000000000000
            1837	0.950000000000000	0	0	0.950000000000000
            1838	0.950000000000000	0	0	0.950000000000000
            1839	0.950000000000000	0	0	0.950000000000000
            1840	0.950000000000000	0	0	0.950000000000000
            1841	0.803796314160256	0	0	0.803796314160256
            1842	0.950000000000000	0	0	0.950000000000000
            1843	0.700720951616440	0	0	0.700720951616440
            1844	0.536779435010483	0	0	0.536779435010483
            1845	0.307933354641913	0	0	0.307933354641913
            1846	0.406333768582621	0	0	0.406333768582621
            1847	0.198463243087802	0	0	0.198463243087802
            1848	0.149413271329718	0	0	0.149413271329718
            1849	0.0606507179611852	0	0	0.0606507179611852
            1850	0.0808722001910162	0	0	0.0808722001910162
            1851	0.0500000000000000	0	0	0.0500000000000000
            1852	0.0500000000000000	0	0	0.0500000000000000
            1853	0.0500000000000000	0	0	0.0500000000000000
            1854	0.0500000000000000	0	0	0.0500000000000000
            1855	0.0500000000000000	0	0	0.0500000000000000
            1856	0.0500000000000000	0	0	0.0500000000000000
            1857	0.0500000000000000	0	0	0.0500000000000000
            1858	0.0500000000000000	0	0	0.0500000000000000
            1859	0.0500000000000000	0	0	0.0500000000000000
            1860	0.0500000000000000	0	0	0.0500000000000000
            1861	0.0500000000000000	0	0	0.0500000000000000
            1862	0.0500000000000000	0	0	0.0500000000000000
            1863	0.0500000000000000	0	0	0.0500000000000000
            1864	0.0500000000000000	0	0	0.0500000000000000
            1865	0.0500000000000000	0	0	0.0500000000000000
            1866	0.0500000000000000	0	0	0.0500000000000000
            1867	0.0500000000000000	0	0	0.0500000000000000
            1868	0.0500000000000000	0	0	0.0500000000000000
            1869	0.0500000000000000	0	0	0.0500000000000000
            1870	0.0600282679539307	0	0	0.0600282679539307
            1871	0.143388481629892	0	0	0.143388481629892
            1872	0.0927451992822336	0	0	0.0927451992822336
            1873	0.180091706102300	0	0	0.180091706102300
            1874	0.282883364374048	0	0	0.282883364374048
            1875	0.459904835102092	0	0	0.459904835102092
            1876	0.289583913940433	0	0	0.289583913940433
            1877	0.384968693678828	0	0	0.384968693678828
            1878	0.615683129018894	0	0	0.615683129018894
            1879	0.678688492752129	0	0	0.678688492752129
            1880	0.423008502060554	0	0	0.423008502060554
            1881	0.385471708578032	0	0	0.385471708578032
            1882	0.617037399965952	0	0	0.617037399965952
            1883	0.467120667150047	0	0	0.467120667150047
            1884	0.296719576855888	0	0	0.296719576855888
            1885	0.210454349877117	0	0	0.210454349877117
            1886	0.309845621285636	0	0	0.309845621285636
            1887	0.222758782912397	0	0	0.222758782912397
            1888	0.182678752084126	0	0	0.182678752084126
            1889	0.253400446353007	0	0	0.253400446353007
            1890	0.250772463464882	0	0	0.250772463464882
            1891	0.398053041333384	0	0	0.398053041333384
            1892	0.433261772859228	0	0	0.433261772859228
            1893	0.683650797910872	0	0	0.683650797910872
            1894	0.620286760113396	0	0	0.620286760113396
            1895	0.825090130670115	0	0	0.825090130670115
            1896	0.909918251962075	0	0	0.909918251962075
            1897	0.950000000000000	0	0	0.950000000000000
            1898	0.910845378925582	0	0	0.910845378925582
            1899	0.835018455532975	0	0	0.835018455532975
            1900	0.913464366927386	0	0	0.913464366927386
            1901	0.694950916680415	0	0	0.694950916680415
            1902	0.649842913079556	0	0	0.649842913079556
            1903	0.463646087857013	0	0	0.463646087857013
            1904	0.459190212748537	0	0	0.459190212748537
            1905	0.296652850568461	0	0	0.296652850568461
            1906	0.362531679810505	0	0	0.362531679810505
            1907	0.359903824254756	0	0	0.359903824254756
            1908	0.228869011438425	0	0	0.228869011438425
            1909	0.219838555789821	0	0	0.219838555789821
            1910	0.403566755281257	0	0	0.403566755281257
            1911	0.422812575632091	0	0	0.422812575632091
            1912	0.217964217066123	0	0	0.217964217066123
            1913	0.192039732325075	0	0	0.192039732325075
            1914	0.378354358382205	0	0	0.378354358382205
            1915	0.281296190509174	0	0	0.281296190509174
            1916	0.142197159653853	0	0	0.142197159653853
            1917	0.0870957605681540	0	0	0.0870957605681540
            1918	0.172473703401559	0	0	0.172473703401559
            1919	0.0870450774901912	0	0	0.0870450774901912
            1920	0.0500000000000000	0	0	0.0500000000000000
            1921	0.458177149616924	0	0	0.458177149616924
            1922	0.613764253312195	0	0	0.613764253312195
            1923	0.822803051791615	0	0	0.822803051791615
            1924	0.614104730956225	0	0	0.614104730956225
            1925	0.677215410078246	0	0	0.677215410078246
            1926	0.907740355910982	0	0	0.907740355910982
            1927	0.825209164661964	0	0	0.825209164661964
            1928	0.614814498794021	0	0	0.614814498794021
            1929	0.460585416665821	0	0	0.460585416665821
            1930	0.621313295750710	0	0	0.621313295750710
            1931	0.395039828692887	0	0	0.395039828692887
            1932	0.287334896373542	0	0	0.287334896373542
            1933	0.154487818612488	0	0	0.154487818612488
            1934	0.227059817946864	0	0	0.227059817946864
            1935	0.140827518681912	0	0	0.140827518681912
            1936	0.0800724300708113	0	0	0.0800724300708113
            1937	0.0503979234160990	0	0	0.0503979234160990
            1938	0.115231177683778	0	0	0.115231177683778
            1939	0.119094999012035	0	0	0.119094999012035
            1940	0.0500000000000000	0	0	0.0500000000000000
            1941	0.0500000000000000	0	0	0.0500000000000000
            1942	0.135115915921437	0	0	0.135115915921437
            1943	0.162305037042778	0	0	0.162305037042778
            1944	0.0590578661516374	0	0	0.0590578661516374
            1945	0.0772465171881632	0	0	0.0772465171881632
            1946	0.203485235075695	0	0	0.203485235075695
            1947	0.252823307059297	0	0	0.252823307059297
            1948	0.102491976517469	0	0	0.102491976517469
            1949	0.135821944336275	0	0	0.135821944336275
            1950	0.297074677454034	0	0	0.297074677454034
            1951	0.331736089144767	0	0	0.331736089144767
            1952	0.188908874379703	0	0	0.188908874379703
            1953	0.289995483285810	0	0	0.289995483285810
            1954	0.376383667627446	0	0	0.376383667627446
            1955	0.468168940124718	0	0	0.468168940124718
            1956	0.471514791950540	0	0	0.471514791950540
            1957	0.739458709324599	0	0	0.739458709324599
            1958	0.631419725896219	0	0	0.631419725896219
            1959	0.850319041005559	0	0	0.850319041005559
            1960	0.950000000000000	0	0	0.950000000000000
            1961	0.950000000000000	0	0	0.950000000000000
            1962	0.950000000000000	0	0	0.950000000000000
            1963	0.950000000000000	0	0	0.950000000000000
            1964	0.950000000000000	0	0	0.950000000000000
            1965	0.950000000000000	0	0	0.950000000000000
            1966	0.950000000000000	0	0	0.950000000000000
            1967	0.950000000000000	0	0	0.950000000000000
            1968	0.950000000000000	0	0	0.950000000000000
            1969	0.950000000000000	0	0	0.950000000000000
            1970	0.950000000000000	0	0	0.950000000000000
            1971	0.672590459615937	0	0	0.672590459615937
            1972	0.755579928392126	0	0	0.755579928392126
            1973	0.442202325648694	0	0	0.442202325648694
            1974	0.396732894002467	0	0	0.396732894002467
            1975	0.195785552165656	0	0	0.195785552165656
            1976	0.217190377116754	0	0	0.217190377116754
            1977	0.0887949851130588	0	0	0.0887949851130588
            1978	0.0802884993796209	0	0	0.0802884993796209
            1979	0.0500000000000000	0	0	0.0500000000000000
            1980	0.0500000000000000	0	0	0.0500000000000000
            1981	0.0500000000000000	0	0	0.0500000000000000
            1982	0.0500000000000000	0	0	0.0500000000000000
            1983	0.0500000000000000	0	0	0.0500000000000000
            1984	0.0500000000000000	0	0	0.0500000000000000
            1985	0.0500000000000000	0	0	0.0500000000000000
            1986	0.0500000000000000	0	0	0.0500000000000000
            1987	0.0500000000000000	0	0	0.0500000000000000
            1988	0.0500000000000000	0	0	0.0500000000000000
            1989	0.0500000000000000	0	0	0.0500000000000000
            1990	0.0500000000000000	0	0	0.0500000000000000
            1991	0.0500000000000000	0	0	0.0500000000000000
            1992	0.0500000000000000	0	0	0.0500000000000000
            1993	0.0500000000000000	0	0	0.0500000000000000
            1994	0.0500000000000000	0	0	0.0500000000000000
            1995	0.0500000000000000	0	0	0.0500000000000000
            1996	0.0500000000000000	0	0	0.0500000000000000
            1997	0.0808437609042273	0	0	0.0808437609042273
            1998	0.0946380105009356	0	0	0.0946380105009356
            1999	0.215302179903391	0	0	0.215302179903391
            2000	0.191664563818850	0	0	0.191664563818850
            2001	0.377741492993200	0	0	0.377741492993200
            2002	0.418758085050846	0	0	0.418758085050846
            2003	0.678446855109046	0	0	0.678446855109046
            2004	0.614467799832816	0	0	0.614467799832816
            2005	0.823149214498381	0	0	0.823149214498381
            2006	0.907768511011669	0	0	0.907768511011669
            2007	0.950000000000000	0	0	0.950000000000000
            2008	0.907591583104925	0	0	0.907591583104925
            2009	0.824204682761207	0	0	0.824204682761207
            2010	0.907879645196934	0	0	0.907879645196934
            2011	0.680169235845261	0	0	0.680169235845261
            2012	0.619661034830040	0	0	0.619661034830040
            2013	0.397084798960402	0	0	0.397084798960402
            2014	0.427641668113074	0	0	0.427641668113074
            2015	0.245813208099992	0	0	0.245813208099992
            2016	0.249380930210861	0	0	0.249380930210861
            2017	0.221111461057642	0	0	0.221111461057642
            2018	0.174136023461474	0	0	0.174136023461474
            2019	0.202241842032844	0	0	0.202241842032844
            2020	0.308214388336237	0	0	0.308214388336237
            2021	0.465809888858422	0	0	0.465809888858422
            2022	0.289672818224701	0	0	0.289672818224701
            2023	0.380296798117004	0	0	0.380296798117004
            2024	0.616814257483670	0	0	0.616814257483670
            2025	0.682650450136875	0	0	0.682650450136875
            2026	0.423184560664974	0	0	0.423184560664974
            2027	0.403625995767649	0	0	0.403625995767649
            2028	0.633859030589430	0	0	0.633859030589430
            2029	0.516736296308964	0	0	0.516736296308964
            2030	0.359912168542161	0	0	0.359912168542161
            2031	0.362533336480927	0	0	0.362533336480927
            2032	0.422795032227834	0	0	0.422795032227834
            2033	0.422793972256521	0	0	0.422793972256521
            2034	0.463644972877215	0	0	0.463644972877215
            2035	0.649826790325716	0	0	0.649826790325716
            2036	0.516727021268540	0	0	0.516727021268540
            2037	0.633796756663943	0	0	0.633796756663943
            2038	0.834928924222353	0	0	0.834928924222353
            2039	0.910436398259006	0	0	0.910436398259006
            2040	0.682286427329574	0	0	0.682286427329574
            2041	0.614987170302124	0	0	0.614987170302124
            2042	0.823385853639776	0	0	0.823385853639776
            2043	0.613837637174354	0	0	0.613837637174354
            2044	0.458082213074038	0	0	0.458082213074038
            2045	0.281002449437902	0	0	0.281002449437902
            2046	0.376630346981410	0	0	0.376630346981410
            2047	0.190113886110251	0	0	0.190113886110251
            2048	0.141835318033456	0	0	0.141835318033456
            2049	0.676717483058458	0	0	0.676717483058458
            2050	0.613833378871192	0	0	0.613833378871192
            2051	0.823380266881494	0	0	0.823380266881494
            2052	0.907372811679945	0	0	0.907372811679945
            2053	0.950000000000000	0	0	0.950000000000000
            2054	0.910432835755944	0	0	0.910432835755944
            2055	0.834947115962319	0	0	0.834947115962319
            2056	0.913209551549088	0	0	0.913209551549088
            2057	0.694961140161731	0	0	0.694961140161731
            2058	0.649979497490744	0	0	0.649979497490744
            2059	0.464482410693799	0	0	0.464482410693799
            2060	0.459603790261374	0	0	0.459603790261374
            2061	0.298553640893940	0	0	0.298553640893940
            2062	0.366247687016444	0	0	0.366247687016444
            2063	0.373554578044723	0	0	0.373554578044723
            2064	0.235901956811315	0	0	0.235901956811315
            2065	0.241551509956877	0	0	0.241551509956877
            2066	0.445360185226728	0	0	0.445360185226728
            2067	0.529968218899325	0	0	0.529968218899325
            2068	0.274248889850346	0	0	0.274248889850346
            2069	0.315243511301851	0	0	0.315243511301851
            2070	0.609681084312658	0	0	0.609681084312658
            2071	0.704166847271492	0	0	0.704166847271492
            2072	0.371218603229176	0	0	0.371218603229176
            2073	0.450433246708492	0	0	0.450433246708492
            2074	0.830448128983105	0	0	0.830448128983105
            2075	0.950000000000000	0	0	0.950000000000000
            2076	0.538026226217868	0	0	0.538026226217868
            2077	0.597550726361194	0	0	0.597550726361194
            2078	0.950000000000000	0	0	0.950000000000000
            2079	0.950000000000000	0	0	0.950000000000000
            2080	0.600746384472998	0	0	0.600746384472998
            2081	0.558082836609093	0	0	0.558082836609093
            2082	0.836975608890190	0	0	0.836975608890190
            2083	0.645056358670307	0	0	0.645056358670307
            2084	0.517842289678639	0	0	0.517842289678639
            2085	0.530599182163974	0	0	0.530599182163974
            2086	0.495228915365785	0	0	0.495228915365785
            2087	0.431346310558983	0	0	0.431346310558983
            2088	0.611715592521764	0	0	0.611715592521764
            2089	0.734094071545232	0	0	0.734094071545232
            2090	0.444397690652926	0	0	0.444397690652926
            2091	0.491990988190720	0	0	0.491990988190720
            2092	0.847186421811194	0	0	0.847186421811194
            2093	0.898189591203893	0	0	0.898189591203893
            2094	0.524241555203013	0	0	0.524241555203013
            2095	0.503621364940016	0	0	0.503621364940016
            2096	0.849872162239219	0	0	0.849872162239219
            2097	0.700883469830582	0	0	0.700883469830582
            2098	0.420849169477922	0	0	0.420849169477922
            2099	0.299414774003589	0	0	0.299414774003589
            2100	0.494030858347319	0	0	0.494030858347319
            2101	0.293388184239599	0	0	0.293388184239599
            2102	0.179063517371535	0	0	0.179063517371535
            2103	0.0895963350451980	0	0	0.0895963350451980
            2104	0.145452979597480	0	0	0.145452979597480
            2105	0.0599539632401141	0	0	0.0599539632401141
            2106	0.0500000000000000	0	0	0.0500000000000000
            2107	0.0500000000000000	0	0	0.0500000000000000
            2108	0.0500000000000000	0	0	0.0500000000000000
            2109	0.0500000000000000	0	0	0.0500000000000000
            2110	0.0500000000000000	0	0	0.0500000000000000
            2111	0.0500000000000000	0	0	0.0500000000000000
            2112	0.0500000000000000	0	0	0.0500000000000000
            2113	0.0500000000000000	0	0	0.0500000000000000
            2114	0.0500000000000000	0	0	0.0500000000000000
            2115	0.0500000000000000	0	0	0.0500000000000000
            2116	0.0500000000000000	0	0	0.0500000000000000
            2117	0.0500000000000000	0	0	0.0500000000000000
            2118	0.0500000000000000	0	0	0.0500000000000000
            2119	0.0603677950887084	0	0	0.0603677950887084
            2120	0.0500000000000000	0	0	0.0500000000000000
            2121	0.0500000000000000	0	0	0.0500000000000000
            2122	0.0846497882655175	0	0	0.0846497882655175
            2123	0.107170581925612	0	0	0.107170581925612
            2124	0.0566970708166817	0	0	0.0566970708166817
            2125	0.105907385739146	0	0	0.105907385739146
            2126	0.137837636703686	0	0	0.137837636703686
            2127	0.200724583253864	0	0	0.200724583253864
            2128	0.210230230959319	0	0	0.210230230959319
            2129	0.388955061334061	0	0	0.388955061334061
            2130	0.317105764481567	0	0	0.317105764481567
            2131	0.476079076102174	0	0	0.476079076102174
            2132	0.619920986998843	0	0	0.619920986998843
            2133	0.825184513646413	0	0	0.825184513646413
            2134	0.621261709762673	0	0	0.621261709762673
            2135	0.679254066491871	0	0	0.679254066491871
            2136	0.907927913876487	0	0	0.907927913876487
            2137	0.823236661587917	0	0	0.823236661587917
            2138	0.614575940158041	0	0	0.614575940158041
            2139	0.458729097902102	0	0	0.458729097902102
            2140	0.615495035508960	0	0	0.615495035508960
            2141	0.382853349686284	0	0	0.382853349686284
            2142	0.283583655945898	0	0	0.283583655945898
            2143	0.149412344876622	0	0	0.149412344876622
            2144	0.208363595495665	0	0	0.208363595495665
            2145	0.122891664929207	0	0	0.122891664929207
            2146	0.0771453942964235	0	0	0.0771453942964235
            2147	0.0563100001747627	0	0	0.0563100001747627
            2148	0.114064022991843	0	0	0.114064022991843
            2149	0.149775925021792	0	0	0.149775925021792
            2150	0.0648782224817825	0	0	0.0648782224817825
            2151	0.0819983008703679	0	0	0.0819983008703679
            2152	0.193792232022504	0	0	0.193792232022504
            2153	0.218298179699765	0	0	0.218298179699765
            2154	0.0948655013314256	0	0	0.0948655013314256
            2155	0.105953366554576	0	0	0.105953366554576
            2156	0.219892466533898	0	0	0.219892466533898
            2157	0.228880279198686	0	0	0.228880279198686
            2158	0.137870421753840	0	0	0.137870421753840
            2159	0.226328265823759	0	0	0.226328265823759
            2160	0.296662507311495	0	0	0.296662507311495
            2161	0.459197657486599	0	0	0.459197657486599
            2162	0.395000120747584	0	0	0.395000120747584
            2163	0.621504254730134	0	0	0.621504254730134
            2164	0.694916916595869	0	0	0.694916916595869
            2165	0.913244251762698	0	0	0.913244251762698
            2166	0.825497997102274	0	0	0.825497997102274
            2167	0.908054630820392	0	0	0.908054630820392
            2168	0.950000000000000	0	0	0.950000000000000
            2169	0.907417314333018	0	0	0.907417314333018
            2170	0.823078438012872	0	0	0.823078438012872
            2171	0.613960516785930	0	0	0.613960516785930
            2172	0.676749605576821	0	0	0.676749605576821
            2173	0.415272689273808	0	0	0.415272689273808
            2174	0.376760836684727	0	0	0.376760836684727
            2175	0.190186842559171	0	0	0.190186842559171
            2176	0.209624964376324	0	0	0.209624964376324
            2177	0.458081792449620	0	0	0.458081792449620
            2178	0.281296210861450	0	0	0.281296210861450
            2179	0.378355014793143	0	0	0.378355014793143
            2180	0.614987062762788	0	0	0.614987062762788
            2181	0.682290764421290	0	0	0.682290764421290
            2182	0.422819274112441	0	0	0.422819274112441
            2183	0.403619514890348	0	0	0.403619514890348
            2184	0.633835632322313	0	0	0.633835632322313
            2185	0.516981513814604	0	0	0.516981513814604
            2186	0.360243900665699	0	0	0.360243900665699
            2187	0.364339134892418	0	0	0.364339134892418
            2188	0.424151727715521	0	0	0.424151727715521
            2189	0.428778109656343	0	0	0.428778109656343
            2190	0.471584493611107	0	0	0.471584493611107
            2191	0.678715471722861	0	0	0.678715471722861
            2192	0.538593115097704	0	0	0.538593115097704
            2193	0.700264896765745	0	0	0.700264896765745
            2194	0.922215441445380	0	0	0.922215441445380
            2195	0.950000000000000	0	0	0.950000000000000
            2196	0.851079677143661	0	0	0.851079677143661
            2197	0.950000000000000	0	0	0.950000000000000
            2198	0.950000000000000	0	0	0.950000000000000
            2199	0.950000000000000	0	0	0.950000000000000
            2200	0.950000000000000	0	0	0.950000000000000
            2201	0.950000000000000	0	0	0.950000000000000
            2202	0.950000000000000	0	0	0.950000000000000
            2203	0.950000000000000	0	0	0.950000000000000
            2204	0.950000000000000	0	0	0.950000000000000
            2205	0.950000000000000	0	0	0.950000000000000
            2206	0.950000000000000	0	0	0.950000000000000
            2207	0.950000000000000	0	0	0.950000000000000
            2208	0.950000000000000	0	0	0.950000000000000
            2209	0.950000000000000	0	0	0.950000000000000
            2210	0.950000000000000	0	0	0.950000000000000
            2211	0.883426215382415	0	0	0.883426215382415
            2212	0.799189739700307	0	0	0.799189739700307
            2213	0.512882223916683	0	0	0.512882223916683
            2214	0.522201547127975	0	0	0.522201547127975
            2215	0.281834347807079	0	0	0.281834347807079
            2216	0.333141266849650	0	0	0.333141266849650
            2217	0.258674381328062	0	0	0.258674381328062
            2218	0.160493199015490	0	0	0.160493199015490
            2219	0.117348015214139	0	0	0.117348015214139
            2220	0.249043563883593	0	0	0.249043563883593
            2221	0.256842549669561	0	0	0.256842549669561
            2222	0.107799288967469	0	0	0.107799288967469
            2223	0.101613867810890	0	0	0.101613867810890
            2224	0.247612375298029	0	0	0.247612375298029
            2225	0.208939690638074	0	0	0.208939690638074
            2226	0.0864295017474654	0	0	0.0864295017474654
            2227	0.0641266947656042	0	0	0.0641266947656042
            2228	0.150185595843129	0	0	0.150185595843129
            2229	0.0913850852257737	0	0	0.0913850852257737
            2230	0.0500000000000000	0	0	0.0500000000000000
            2231	0.0500000000000000	0	0	0.0500000000000000
            2232	0.0500000000000000	0	0	0.0500000000000000
            2233	0.0500000000000000	0	0	0.0500000000000000
            2234	0.0500000000000000	0	0	0.0500000000000000
            2235	0.0500000000000000	0	0	0.0500000000000000
            2236	0.0500000000000000	0	0	0.0500000000000000
            2237	0.0500000000000000	0	0	0.0500000000000000
            2238	0.0500000000000000	0	0	0.0500000000000000
            2239	0.0500000000000000	0	0	0.0500000000000000
            2240	0.0500000000000000	0	0	0.0500000000000000
            2241	0.0500000000000000	0	0	0.0500000000000000
            2242	0.0562931012848196	0	0	0.0562931012848196
            2243	0.0994047347814918	0	0	0.0994047347814918
            2244	0.0500000000000000	0	0	0.0500000000000000
            2245	0.0889440089431050	0	0	0.0889440089431050
            2246	0.178738882395619	0	0	0.178738882395619
            2247	0.283967307263435	0	0	0.283967307263435
            2248	0.143360175259821	0	0	0.143360175259821
            2249	0.193827504641733	0	0	0.193827504641733
            2250	0.379252221513302	0	0	0.379252221513302
            2251	0.421730819107159	0	0	0.421730819107159
            2252	0.222027327590113	0	0	0.222027327590113
            2253	0.226268199331111	0	0	0.226268199331111
            2254	0.394916769484220	0	0	0.394916769484220
            2255	0.324919201735957	0	0	0.324919201735957
            2256	0.228872493321645	0	0	0.228872493321645
            2257	0.259461044492678	0	0	0.259461044492678
            2258	0.259466321755500	0	0	0.259466321755500
            2259	0.228871002921380	0	0	0.228871002921380
            2260	0.324900808964338	0	0	0.324900808964338
            2261	0.394849555199691	0	0	0.394849555199691
            2262	0.226252178439661	0	0	0.226252178439661
            2263	0.221959405626064	0	0	0.221959405626064
            2264	0.421479412090652	0	0	0.421479412090652
            2265	0.378411752209366	0	0	0.378411752209366
            2266	0.193593876038039	0	0	0.193593876038039
            2267	0.142701547664597	0	0	0.142701547664597
            2268	0.281616325734749	0	0	0.281616325734749
            2269	0.173382491049847	0	0	0.173382491049847
            2270	0.0874391179451545	0	0	0.0874391179451545
            2271	0.0500000000000000	0	0	0.0500000000000000
            2272	0.0896385800282997	0	0	0.0896385800282997
            2273	0.0500000000000000	0	0	0.0500000000000000
            2274	0.0500000000000000	0	0	0.0500000000000000
            2275	0.0500000000000000	0	0	0.0500000000000000
            2276	0.0500000000000000	0	0	0.0500000000000000
            2277	0.0500000000000000	0	0	0.0500000000000000
            2278	0.0500000000000000	0	0	0.0500000000000000
            2279	0.0500000000000000	0	0	0.0500000000000000
            2280	0.0500000000000000	0	0	0.0500000000000000
            2281	0.0500000000000000	0	0	0.0500000000000000
            2282	0.0500000000000000	0	0	0.0500000000000000
            2283	0.0500000000000000	0	0	0.0500000000000000
            2284	0.0500000000000000	0	0	0.0500000000000000
            2285	0.0791799891039228	0	0	0.0791799891039228
            2286	0.0500000000000000	0	0	0.0500000000000000
            2287	0.0922577751553427	0	0	0.0922577751553427
            2288	0.154579794281411	0	0	0.154579794281411
            2289	0.287996498655901	0	0	0.287996498655901
            2290	0.177640881902862	0	0	0.177640881902862
            2291	0.287342108901817	0	0	0.287342108901817
            2292	0.461727156872248	0	0	0.461727156872248
            2293	0.616314923504784	0	0	0.616314923504784
            2294	0.384427450981785	0	0	0.384427450981785
            2295	0.423659416798763	0	0	0.423659416798763
            2296	0.678812623056731	0	0	0.678812623056731
            2297	0.615488487965833	0	0	0.615488487965833
            2298	0.384195252904429	0	0	0.384195252904429
            2299	0.286618157651453	0	0	0.286618157651453
            2300	0.459150482919055	0	0	0.459150482919055
            2301	0.281766594297478	0	0	0.281766594297478
            2302	0.175890501147326	0	0	0.175890501147326
            2303	0.0887892751990465	0	0	0.0887892751990465
            2304	0.142235023729481	0	0	0.142235023729481
            2305	0.142197224829607	0	0	0.142197224829607
            2306	0.0592260193782678	0	0	0.0592260193782678
            2307	0.0806919293638508	0	0	0.0806919293638508
            2308	0.192040513192593	0	0	0.192040513192593
            2309	0.217971650032665	0	0	0.217971650032665
            2310	0.0946207866994392	0	0	0.0946207866994392
            2311	0.105956180008340	0	0	0.105956180008340
            2312	0.219896589145203	0	0	0.219896589145203
            2313	0.229242219676711	0	0	0.229242219676711
            2314	0.138171324431348	0	0	0.138171324431348
            2315	0.228040487462702	0	0	0.228040487462702
            2316	0.298632260280384	0	0	0.298632260280384
            2317	0.467859432924866	0	0	0.467859432924866
            2318	0.402671683894308	0	0	0.402671683894308
            2319	0.649525367844246	0	0	0.649525367844246
            2320	0.726358141980085	0	0	0.726358141980085
            2321	0.950000000000000	0	0	0.950000000000000
            2322	0.909757317804689	0	0	0.909757317804689
            2323	0.950000000000000	0	0	0.950000000000000
            2324	0.950000000000000	0	0	0.950000000000000
            2325	0.950000000000000	0	0	0.950000000000000
            2326	0.950000000000000	0	0	0.950000000000000
            2327	0.950000000000000	0	0	0.950000000000000
            2328	0.950000000000000	0	0	0.950000000000000
            2329	0.950000000000000	0	0	0.950000000000000
            2330	0.950000000000000	0	0	0.950000000000000
            2331	0.950000000000000	0	0	0.950000000000000
            2332	0.950000000000000	0	0	0.950000000000000
            2333	0.950000000000000	0	0	0.950000000000000
            2334	0.950000000000000	0	0	0.950000000000000
            2335	0.950000000000000	0	0	0.950000000000000
            2336	0.950000000000000	0	0	0.950000000000000
            2337	0.950000000000000	0	0	0.950000000000000
            2338	0.950000000000000	0	0	0.950000000000000
            2339	0.650743413532939	0	0	0.650743413532939
            2340	0.829574322642020	0	0	0.829574322642020
            2341	0.472248265778516	0	0	0.472248265778516
            2342	0.362269906240049	0	0	0.362269906240049
            2343	0.172348678022476	0	0	0.172348678022476
            2344	0.233703782767140	0	0	0.233703782767140
            2345	0.108003051927914	0	0	0.108003051927914
            2346	0.0718863602963426	0	0	0.0718863602963426
            2347	0.0500000000000000	0	0	0.0500000000000000
            2348	0.0566178780744056	0	0	0.0566178780744056
            2349	0.0500000000000000	0	0	0.0500000000000000
            2350	0.0500000000000000	0	0	0.0500000000000000
            2351	0.0500000000000000	0	0	0.0500000000000000
            2352	0.0500000000000000	0	0	0.0500000000000000
            2353	0.0500000000000000	0	0	0.0500000000000000
            2354	0.0500000000000000	0	0	0.0500000000000000
            2355	0.0500000000000000	0	0	0.0500000000000000
            2356	0.0500000000000000	0	0	0.0500000000000000
            2357	0.0500000000000000	0	0	0.0500000000000000
            2358	0.0500000000000000	0	0	0.0500000000000000
            2359	0.0888587208140561	0	0	0.0888587208140561
            2360	0.0500000000000000	0	0	0.0500000000000000
            2361	0.0617454929859623	0	0	0.0617454929859623
            2362	0.142896190342442	0	0	0.142896190342442
            2363	0.192115993149650	0	0	0.192115993149650
            2364	0.0811343934658851	0	0	0.0811343934658851
            2365	0.0929416897940644	0	0	0.0929416897940644
            2366	0.217275918390841	0	0	0.217275918390841
            2367	0.217138212851106	0	0	0.217138212851106
            2368	0.0991210961593004	0	0	0.0991210961593004
            2369	0.117815576705928	0	0	0.117815576705928
            2370	0.220890069321711	0	0	0.220890069321711
            2371	0.277385839706659	0	0	0.277385839706659
            2372	0.178014965877482	0	0	0.178014965877482
            2373	0.299302879777832	0	0	0.299302879777832
            2374	0.420957065275041	0	0	0.420957065275041
            2375	0.632566121841270	0	0	0.632566121841270
            2376	0.465668294283727	0	0	0.465668294283727
            2377	0.617104457303832	0	0	0.617104457303832
            2378	0.829758864716845	0	0	0.829758864716845
            2379	0.910295781940214	0	0	0.910295781940214
            2380	0.680057954678451	0	0	0.680057954678451
            2381	0.621514859475030	0	0	0.621514859475030
            2382	0.826151184305064	0	0	0.826151184305064
            2383	0.620389535499150	0	0	0.620389535499150
            2384	0.476172207872976	0	0	0.476172207872976
            2385	0.317145687302476	0	0	0.317145687302476
            2386	0.389189478144513	0	0	0.389189478144513
            2387	0.210335489750708	0	0	0.210335489750708
            2388	0.200740946942346	0	0	0.200740946942346
            2389	0.137840357687931	0	0	0.137840357687931
            2390	0.105946638791335	0	0	0.105946638791335
            2391	0.0567065698591016	0	0	0.0567065698591016
            2392	0.107157061898628	0	0	0.107157061898628
            2393	0.0845964889811608	0	0	0.0845964889811608
            2394	0.0500000000000000	0	0	0.0500000000000000
            2395	0.0500000000000000	0	0	0.0500000000000000
            2396	0.0602155741087644	0	0	0.0602155741087644
            2397	0.0500000000000000	0	0	0.0500000000000000
            2398	0.0500000000000000	0	0	0.0500000000000000
            2399	0.0500000000000000	0	0	0.0500000000000000
            2400	0.0500000000000000	0	0	0.0500000000000000
            2401	0.0500000000000000	0	0	0.0500000000000000
            2402	0.0500000000000000	0	0	0.0500000000000000
            2403	0.0500000000000000	0	0	0.0500000000000000
            2404	0.0500000000000000	0	0	0.0500000000000000
            2405	0.0500000000000000	0	0	0.0500000000000000
            2406	0.0500000000000000	0	0	0.0500000000000000
            2407	0.0500000000000000	0	0	0.0500000000000000
            2408	0.0500000000000000	0	0	0.0500000000000000
            2409	0.0500000000000000	0	0	0.0500000000000000
            2410	0.0500000000000000	0	0	0.0500000000000000
            2411	0.0500000000000000	0	0	0.0500000000000000
            2412	0.0500000000000000	0	0	0.0500000000000000
            2413	0.0500000000000000	0	0	0.0500000000000000
            2414	0.0500000000000000	0	0	0.0500000000000000
            2415	0.0500000000000000	0	0	0.0500000000000000
            2416	0.0509701796140675	0	0	0.0509701796140675
            2417	0.0997876843959027	0	0	0.0997876843959027
            2418	0.0723600023618370	0	0	0.0723600023618370
            2419	0.117818157223499	0	0	0.117818157223499
            2420	0.162113977312256	0	0	0.162113977312256
            2421	0.217129892848347	0	0	0.217129892848347
            2422	0.157891658387016	0	0	0.157891658387016
            2423	0.174079805647178	0	0	0.174079805647178
            2424	0.239356728563456	0	0	0.239356728563456
            2425	0.217076228587210	0	0	0.217076228587210
            2426	0.157881456309689	0	0	0.157881456309689
            2427	0.117786348925668	0	0	0.117786348925668
            2428	0.161946661782029	0	0	0.161946661782029
            2429	0.0993831440982693	0	0	0.0993831440982693
            2430	0.0722830946037744	0	0	0.0722830946037744
            2431	0.0500000000000000	0	0	0.0500000000000000
            2432	0.0501685302208489	0	0	0.0501685302208489
            2433	0.0500000000000000	0	0	0.0500000000000000
            2434	0.0500000000000000	0	0	0.0500000000000000
            2435	0.0500000000000000	0	0	0.0500000000000000
            2436	0.0500000000000000	0	0	0.0500000000000000
            2437	0.0500000000000000	0	0	0.0500000000000000
            2438	0.0500000000000000	0	0	0.0500000000000000
            2439	0.0500000000000000	0	0	0.0500000000000000
            2440	0.0500000000000000	0	0	0.0500000000000000
            2441	0.0792602877078656	0	0	0.0792602877078656
            2442	0.0500000000000000	0	0	0.0500000000000000
            2443	0.0913203909654090	0	0	0.0913203909654090
            2444	0.155504866036361	0	0	0.155504866036361
            2445	0.292999065233628	0	0	0.292999065233628
            2446	0.177717794441230	0	0	0.177717794441230
            2447	0.294439618852630	0	0	0.294439618852630
            2448	0.481338333942728	0	0	0.481338333942728
            2449	0.676925581135249	0	0	0.676925581135249
            2450	0.414744155140072	0	0	0.414744155140072
            2451	0.508576548625854	0	0	0.508576548625854
            2452	0.831056996247931	0	0	0.831056996247931
            2453	0.930842491498257	0	0	0.930842491498257
            2454	0.567285435635635	0	0	0.567285435635635
            2455	0.605866605351306	0	0	0.605866605351306
            2456	0.950000000000000	0	0	0.950000000000000
            2457	0.950000000000000	0	0	0.950000000000000
            2458	0.636064035136194	0	0	0.636064035136194
            2459	0.644105917527844	0	0	0.644105917527844
            2460	0.950000000000000	0	0	0.950000000000000
            2461	0.950000000000000	0	0	0.950000000000000
            2462	0.602856107084080	0	0	0.602856107084080
            2463	0.502925049583847	0	0	0.502925049583847
            2464	0.893221656112915	0	0	0.893221656112915
            2465	0.663311596587327	0	0	0.663311596587327
            2466	0.366014393551961	0	0	0.366014393551961
            2467	0.229704909754622	0	0	0.229704909754622
            2468	0.424309472048987	0	0	0.424309472048987
            2469	0.231979722801860	0	0	0.231979722801860
            2470	0.123403598194186	0	0	0.123403598194186
            2471	0.0564394032747281	0	0	0.0564394032747281
            2472	0.107933556972086	0	0	0.107933556972086
            2473	0.0500000000000000	0	0	0.0500000000000000
            2474	0.0500000000000000	0	0	0.0500000000000000
            2475	0.0500000000000000	0	0	0.0500000000000000
            2476	0.0500000000000000	0	0	0.0500000000000000
            2477	0.0500000000000000	0	0	0.0500000000000000
            2478	0.0500000000000000	0	0	0.0500000000000000
            2479	0.0500000000000000	0	0	0.0500000000000000
            2480	0.0500000000000000	0	0	0.0500000000000000
            2481	0.0500000000000000	0	0	0.0500000000000000
            2482	0.0500000000000000	0	0	0.0500000000000000
            2483	0.0591626197025610	0	0	0.0591626197025610
            2484	0.0500000000000000	0	0	0.0500000000000000
            2485	0.0878631989850946	0	0	0.0878631989850946
            2486	0.142003978879964	0	0	0.142003978879964
            2487	0.281128403346589	0	0	0.281128403346589
            2488	0.172892867540699	0	0	0.172892867540699
            2489	0.281515515077808	0	0	0.281515515077808
            2490	0.458294398804784	0	0	0.458294398804784
            2491	0.615680314791067	0	0	0.615680314791067
            2492	0.378635114692509	0	0	0.378635114692509
            2493	0.423718141923662	0	0	0.423718141923662
            2494	0.684779007629181	0	0	0.684779007629181
            2495	0.641952380543463	0	0	0.641952380543463
            2496	0.406613868732978	0	0	0.406613868732978
            2497	0.368632490173854	0	0	0.368632490173854
            2498	0.539618741263564	0	0	0.539618741263564
            2499	0.476624793308945	0	0	0.476624793308945
            2500	0.383291893335986	0	0	0.383291893335986
            2501	0.504424914932129	0	0	0.504424914932129
            2502	0.527802841321928	0	0	0.527802841321928
            2503	0.685956951541977	0	0	0.685956951541977
            2504	0.715872395456880	0	0	0.715872395456880
            2505	0.923068382923923	0	0	0.923068382923923
            2506	0.858714903945949	0	0	0.858714903945949
            2507	0.928630247301521	0	0	0.928630247301521
            2508	0.950000000000000	0	0	0.950000000000000
            2509	0.911547991900772	0	0	0.911547991900772
            2510	0.837290184288364	0	0	0.837290184288364
            2511	0.623456482172029	0	0	0.623456482172029
            2512	0.680645368446238	0	0	0.680645368446238
            2513	0.419935409024100	0	0	0.419935409024100
            2514	0.382663831208582	0	0	0.382663831208582
            2515	0.193898475330879	0	0	0.193898475330879
            2516	0.215838942566288	0	0	0.215838942566288
            2517	0.0948397905167616	0	0	0.0948397905167616
            2518	0.0816671060626292	0	0	0.0816671060626292
            2519	0.0500000000000000	0	0	0.0500000000000000
            2520	0.0500000000000000	0	0	0.0500000000000000
            2521	0.0500000000000000	0	0	0.0500000000000000
            2522	0.0500000000000000	0	0	0.0500000000000000
            2523	0.0500000000000000	0	0	0.0500000000000000
            2524	0.0500000000000000	0	0	0.0500000000000000
            2525	0.0500000000000000	0	0	0.0500000000000000
            2526	0.0500000000000000	0	0	0.0500000000000000
            2527	0.0500000000000000	0	0	0.0500000000000000
            2528	0.0500000000000000	0	0	0.0500000000000000
            2529	0.0500000000000000	0	0	0.0500000000000000
            2530	0.0500000000000000	0	0	0.0500000000000000
            2531	0.0500000000000000	0	0	0.0500000000000000
            2532	0.0500000000000000	0	0	0.0500000000000000
            2533	0.0500000000000000	0	0	0.0500000000000000
            2534	0.0500000000000000	0	0	0.0500000000000000
            2535	0.0500000000000000	0	0	0.0500000000000000
            2536	0.0500000000000000	0	0	0.0500000000000000
            2537	0.0500000000000000	0	0	0.0500000000000000
            2538	0.0500000000000000	0	0	0.0500000000000000
            2539	0.0500000000000000	0	0	0.0500000000000000
            2540	0.0500000000000000	0	0	0.0500000000000000
            2541	0.0500000000000000	0	0	0.0500000000000000
            2542	0.0500000000000000	0	0	0.0500000000000000
            2543	0.0887924442392852	0	0	0.0887924442392852
            2544	0.0501923674963886	0	0	0.0501923674963886
            2545	0.0993950917725213	0	0	0.0993950917725213
            2546	0.175891893702919	0	0	0.175891893702919
            2547	0.286616263297945	0	0	0.286616263297945
            2548	0.161950994371099	0	0	0.161950994371099
            2549	0.217074269526866	0	0	0.217074269526866
            2550	0.384180056105047	0	0	0.384180056105047
            2551	0.423590490018140	0	0	0.423590490018140
            2552	0.239341150908358	0	0	0.239341150908358
            2553	0.217072674129397	0	0	0.217072674129397
            2554	0.384179850882782	0	0	0.384179850882782
            2555	0.286615623443772	0	0	0.286615623443772
            2556	0.161946020184339	0	0	0.161946020184339
            2557	0.0993830647153189	0	0	0.0993830647153189
            2558	0.175890346553041	0	0	0.175890346553041
            2559	0.0887893773252589	0	0	0.0887893773252589
            2560	0.0501685305850643	0	0	0.0501685305850643
            2561	0.0500000000000000	0	0	0.0500000000000000
            2562	0.0500000000000000	0	0	0.0500000000000000
            2563	0.0500000000000000	0	0	0.0500000000000000
            2564	0.0500000000000000	0	0	0.0500000000000000
            2565	0.0500000000000000	0	0	0.0500000000000000
            2566	0.0500000000000000	0	0	0.0500000000000000
            2567	0.0500000000000000	0	0	0.0500000000000000
            2568	0.0500000000000000	0	0	0.0500000000000000
            2569	0.0500000000000000	0	0	0.0500000000000000
            2570	0.0500000000000000	0	0	0.0500000000000000
            2571	0.0500000000000000	0	0	0.0500000000000000
            2572	0.0500000000000000	0	0	0.0500000000000000
            2573	0.0892258015421630	0	0	0.0892258015421630
            2574	0.0500000000000000	0	0	0.0500000000000000
            2575	0.0615700216905282	0	0	0.0615700216905282
            2576	0.148397788528122	0	0	0.148397788528122
            2577	0.209137912772684	0	0	0.209137912772684
            2578	0.0867752825038501	0	0	0.0867752825038501
            2579	0.106169492530965	0	0	0.106169492530965
            2580	0.256152028033447	0	0	0.256152028033447
            2581	0.284727045397004	0	0	0.284727045397004
            2582	0.117670198653165	0	0	0.117670198653165
            2583	0.124017928491196	0	0	0.124017928491196
            2584	0.301905262392950	0	0	0.301905262392950
            2585	0.313289302103916	0	0	0.313289302103916
            2586	0.127426595971068	0	0	0.127426595971068
            2587	0.125398404699724	0	0	0.125398404699724
            2588	0.312441505472890	0	0	0.312441505472890
            2589	0.287294119971140	0	0	0.287294119971140
            2590	0.113479304236963	0	0	0.113479304236963
            2591	0.0912099872547312	0	0	0.0912099872547312
            2592	0.235100862784284	0	0	0.235100862784284
            2593	0.167731616501807	0	0	0.167731616501807
            2594	0.0638320105293979	0	0	0.0638320105293979
            2595	0.0500000000000000	0	0	0.0500000000000000
            2596	0.103239517589842	0	0	0.103239517589842
            2597	0.0544733944470028	0	0	0.0544733944470028
            2598	0.0500000000000000	0	0	0.0500000000000000
            2599	0.0500000000000000	0	0	0.0500000000000000
            2600	0.0500000000000000	0	0	0.0500000000000000
            2601	0.0500000000000000	0	0	0.0500000000000000
            2602	0.0500000000000000	0	0	0.0500000000000000
            2603	0.0500000000000000	0	0	0.0500000000000000
            2604	0.0500000000000000	0	0	0.0500000000000000
            2605	0.0500000000000000	0	0	0.0500000000000000
            2606	0.0500000000000000	0	0	0.0500000000000000
            2607	0.0500000000000000	0	0	0.0500000000000000
            2608	0.0500000000000000	0	0	0.0500000000000000
            2609	0.0500000000000000	0	0	0.0500000000000000
            2610	0.0500000000000000	0	0	0.0500000000000000
            2611	0.0871160095395409	0	0	0.0871160095395409
            2612	0.0790019018157894	0	0	0.0790019018157894
            2613	0.190198292058611	0	0	0.190198292058611
            2614	0.209871727229409	0	0	0.209871727229409
            2615	0.416026190872902	0	0	0.416026190872902
            2616	0.376835084091755	0	0	0.376835084091755
            2617	0.614406535985955	0	0	0.614406535985955
            2618	0.678737695354069	0	0	0.678737695354069
            2619	0.912125870601273	0	0	0.912125870601273
            2620	0.824992791622105	0	0	0.824992791622105
            2621	0.915051257120883	0	0	0.915051257120883
            2622	0.950000000000000	0	0	0.950000000000000
            2623	0.937044235080418	0	0	0.937044235080418
            2624	0.847915722208021	0	0	0.847915722208021
            2625	0.683322999335374	0	0	0.683322999335374
            2626	0.746948113999563	0	0	0.746948113999563
            2627	0.564236646241379	0	0	0.564236646241379
            2628	0.538760727849365	0	0	0.538760727849365
            2629	0.504839370943642	0	0	0.504839370943642
            2630	0.482719473115062	0	0	0.482719473115062
            2631	0.509303645967635	0	0	0.509303645967635
            2632	0.584655060142714	0	0	0.584655060142714
            2633	0.697743635147540	0	0	0.697743635147540
            2634	0.574062731719512	0	0	0.574062731719512
            2635	0.590297108445515	0	0	0.590297108445515
            2636	0.741080902575091	0	0	0.741080902575091
            2637	0.662599453295649	0	0	0.662599453295649
            2638	0.516620699990405	0	0	0.516620699990405
            2639	0.376446730288381	0	0	0.376446730288381
            2640	0.490719724901232	0	0	0.490719724901232
            2641	0.299747708974807	0	0	0.299747708974807
            2642	0.226904739075814	0	0	0.226904739075814
            2643	0.112907558626693	0	0	0.112907558626693
            2644	0.150905050650419	0	0	0.150905050650419
            2645	0.0626981855912337	0	0	0.0626981855912337
            2646	0.0500000000000000	0	0	0.0500000000000000
            2647	0.0500000000000000	0	0	0.0500000000000000
            2648	0.0500000000000000	0	0	0.0500000000000000
            2649	0.0500000000000000	0	0	0.0500000000000000
            2650	0.0500000000000000	0	0	0.0500000000000000
            2651	0.0500000000000000	0	0	0.0500000000000000
            2652	0.0500000000000000	0	0	0.0500000000000000
            2653	0.0500000000000000	0	0	0.0500000000000000
            2654	0.0500000000000000	0	0	0.0500000000000000
            2655	0.0500000000000000	0	0	0.0500000000000000
            2656	0.0500000000000000	0	0	0.0500000000000000
            2657	0.0500000000000000	0	0	0.0500000000000000
            2658	0.0500000000000000	0	0	0.0500000000000000
            2659	0.0500000000000000	0	0	0.0500000000000000
            2660	0.0500000000000000	0	0	0.0500000000000000
            2661	0.0500000000000000	0	0	0.0500000000000000
            2662	0.0500000000000000	0	0	0.0500000000000000
            2663	0.0500000000000000	0	0	0.0500000000000000
            2664	0.0500000000000000	0	0	0.0500000000000000
            2665	0.0500000000000000	0	0	0.0500000000000000
            2666	0.0500000000000000	0	0	0.0500000000000000
            2667	0.0500000000000000	0	0	0.0500000000000000
            2668	0.0500000000000000	0	0	0.0500000000000000
            2669	0.0590623028643707	0	0	0.0590623028643707
            2670	0.0789773561867393	0	0	0.0789773561867393
            2671	0.190198258064899	0	0	0.190198258064899
            2672	0.142236676347905	0	0	0.142236676347905
            2673	0.281768117335572	0	0	0.281768117335572
            2674	0.376779686541975	0	0	0.376779686541975
            2675	0.613967412482100	0	0	0.613967412482100
            2676	0.459144635580067	0	0	0.459144635580067
            2677	0.615437815882850	0	0	0.615437815882850
            2678	0.822962467733103	0	0	0.822962467733103
            2679	0.907385006577980	0	0	0.907385006577980
            2680	0.678571689266959	0	0	0.678571689266959
            2681	0.615437794168827	0	0	0.615437794168827
            2682	0.822962465843717	0	0	0.822962465843717
            2683	0.613967406590008	0	0	0.613967406590008
            2684	0.459144567876586	0	0	0.459144567876586
            2685	0.281767953619129	0	0	0.281767953619129
            2686	0.376779672291337	0	0	0.376779672291337
            2687	0.190198229788806	0	0	0.190198229788806
            2688	0.142236351705875	0	0	0.142236351705875
            2689	0.0500000000000000	0	0	0.0500000000000000
            2690	0.0500000000000000	0	0	0.0500000000000000
            2691	0.0500000000000000	0	0	0.0500000000000000
            2692	0.0500000000000000	0	0	0.0500000000000000
            2693	0.0500000000000000	0	0	0.0500000000000000
            2694	0.0500000000000000	0	0	0.0500000000000000
            2695	0.0500000000000000	0	0	0.0500000000000000
            2696	0.0500000000000000	0	0	0.0500000000000000
            2697	0.0500000000000000	0	0	0.0500000000000000
            2698	0.0500000000000000	0	0	0.0500000000000000
            2699	0.0500000000000000	0	0	0.0500000000000000
            2700	0.0500000000000000	0	0	0.0500000000000000
            2701	0.0500000000000000	0	0	0.0500000000000000
            2702	0.0500000000000000	0	0	0.0500000000000000
            2703	0.0500000000000000	0	0	0.0500000000000000
            2704	0.0500000000000000	0	0	0.0500000000000000
            2705	0.0500000000000000	0	0	0.0500000000000000
            2706	0.0500000000000000	0	0	0.0500000000000000
            2707	0.0500000000000000	0	0	0.0500000000000000
            2708	0.0500000000000000	0	0	0.0500000000000000
            2709	0.0500000000000000	0	0	0.0500000000000000
            2710	0.0500000000000000	0	0	0.0500000000000000
            2711	0.0500000000000000	0	0	0.0500000000000000
            2712	0.0500000000000000	0	0	0.0500000000000000
            2713	0.0500000000000000	0	0	0.0500000000000000
            2714	0.0500000000000000	0	0	0.0500000000000000
            2715	0.0500000000000000	0	0	0.0500000000000000
            2716	0.0500000000000000	0	0	0.0500000000000000
            2717	0.0500000000000000	0	0	0.0500000000000000
            2718	0.0500000000000000	0	0	0.0500000000000000
            2719	0.0500000000000000	0	0	0.0500000000000000
            2720	0.0500000000000000	0	0	0.0500000000000000
            2721	0.0500000000000000	0	0	0.0500000000000000
            2722	0.0500000000000000	0	0	0.0500000000000000
            2723	0.0500000000000000	0	0	0.0500000000000000
            2724	0.0500000000000000	0	0	0.0500000000000000
            2725	0.0500000000000000	0	0	0.0500000000000000
            2726	0.0500000000000000	0	0	0.0500000000000000
            2727	0.0500000000000000	0	0	0.0500000000000000
            2728	0.0500000000000000	0	0	0.0500000000000000
            2729	0.0500000000000000	0	0	0.0500000000000000
            2730	0.0500000000000000	0	0	0.0500000000000000
            2731	0.0500000000000000	0	0	0.0500000000000000
            2732	0.0500000000000000	0	0	0.0500000000000000
            2733	0.0500000000000000	0	0	0.0500000000000000
            2734	0.0500000000000000	0	0	0.0500000000000000
            2735	0.0500000000000000	0	0	0.0500000000000000
            2736	0.0500000000000000	0	0	0.0500000000000000
            2737	0.0500000000000000	0	0	0.0500000000000000
            2738	0.0500000000000000	0	0	0.0500000000000000
            2739	0.0596237504600327	0	0	0.0596237504600327
            2740	0.0791916069025149	0	0	0.0791916069025149
            2741	0.191008614801107	0	0	0.191008614801107
            2742	0.144445332876736	0	0	0.144445332876736
            2743	0.288720211828869	0	0	0.288720211828869
            2744	0.379292191066172	0	0	0.379292191066172
            2745	0.620470401403840	0	0	0.620470401403840
            2746	0.477063537795016	0	0	0.477063537795016
            2747	0.654259811418287	0	0	0.654259811418287
            2748	0.837482387747797	0	0	0.837482387747797
            2749	0.936809692664000	0	0	0.936809692664000
            2750	0.752067034343450	0	0	0.752067034343450
            2751	0.743694408246088	0	0	0.743694408246088
            2752	0.880248190013636	0	0	0.880248190013636
            2753	0.723834867544584	0	0	0.723834867544584
            2754	0.675627310891479	0	0	0.675627310891479
            2755	0.637607888936896	0	0	0.637607888936896
            2756	0.578478914559475	0	0	0.578478914559475
            2757	0.526131964355550	0	0	0.526131964355550
            2758	0.690249717012425	0	0	0.690249717012425
            2759	0.811494491395979	0	0	0.811494491395979
            2760	0.565054642880116	0	0	0.565054642880116
            2761	0.623024805983487	0	0	0.623024805983487
            2762	0.909140684645040	0	0	0.909140684645040
            2763	0.892106147466799	0	0	0.892106147466799
            2764	0.619495526900725	0	0	0.619495526900725
            2765	0.524684020253762	0	0	0.524684020253762
            2766	0.740400976077666	0	0	0.740400976077666
            2767	0.513253524570459	0	0	0.513253524570459
            2768	0.371281408104833	0	0	0.371281408104833
            2769	0.218126176708645	0	0	0.218126176708645
            2770	0.295808709888662	0	0	0.295808709888662
            2771	0.141461608915025	0	0	0.141461608915025
            2772	0.106140701983904	0	0	0.106140701983904
            2773	0.0500000000000000	0	0	0.0500000000000000
            2774	0.0560681252196355	0	0	0.0560681252196355
            2775	0.0500000000000000	0	0	0.0500000000000000
            2776	0.0500000000000000	0	0	0.0500000000000000
            2777	0.0500000000000000	0	0	0.0500000000000000
            2778	0.0500000000000000	0	0	0.0500000000000000
            2779	0.0500000000000000	0	0	0.0500000000000000
            2780	0.0500000000000000	0	0	0.0500000000000000
            2781	0.0500000000000000	0	0	0.0500000000000000
            2782	0.0500000000000000	0	0	0.0500000000000000
            2783	0.0500000000000000	0	0	0.0500000000000000
            2784	0.0500000000000000	0	0	0.0500000000000000
            2785	0.0500000000000000	0	0	0.0500000000000000
            2786	0.0500000000000000	0	0	0.0500000000000000
            2787	0.0500000000000000	0	0	0.0500000000000000
            2788	0.0500000000000000	0	0	0.0500000000000000
            2789	0.0500000000000000	0	0	0.0500000000000000
            2790	0.0500000000000000	0	0	0.0500000000000000
            2791	0.0500000000000000	0	0	0.0500000000000000
            2792	0.0500000000000000	0	0	0.0500000000000000
            2793	0.0500000000000000	0	0	0.0500000000000000
            2794	0.0500000000000000	0	0	0.0500000000000000
            2795	0.0500000000000000	0	0	0.0500000000000000
            2796	0.0500000000000000	0	0	0.0500000000000000
            2797	0.0870753090625809	0	0	0.0870753090625809
            2798	0.0791090381725898	0	0	0.0791090381725898
            2799	0.190515465963360	0	0	0.190515465963360
            2800	0.209700345761103	0	0	0.209700345761103
            2801	0.415413051278483	0	0	0.415413051278483
            2802	0.377408111888569	0	0	0.377408111888569
            2803	0.614991457039935	0	0	0.614991457039935
            2804	0.676921001950818	0	0	0.676921001950818
            2805	0.907345521344519	0	0	0.907345521344519
            2806	0.824335103867393	0	0	0.824335103867393
            2807	0.908898454468032	0	0	0.908898454468032
            2808	0.950000000000000	0	0	0.950000000000000
            2809	0.907345521209195	0	0	0.907345521209195
            2810	0.824335103858806	0	0	0.824335103858806
            2811	0.614991457013707	0	0	0.614991457013707
            2812	0.676921001528824	0	0	0.676921001528824
            2813	0.415413050257561	0	0	0.415413050257561
            2814	0.377408111823682	0	0	0.377408111823682
            2815	0.190515465817909	0	0	0.190515465817909
            2816	0.209700343731927	0	0	0.209700343731927
            2817	0.0500000000000000	0	0	0.0500000000000000
            2818	0.0500000000000000	0	0	0.0500000000000000
            2819	0.0500000000000000	0	0	0.0500000000000000
            2820	0.0500000000000000	0	0	0.0500000000000000
            2821	0.0500000000000000	0	0	0.0500000000000000
            2822	0.0500000000000000	0	0	0.0500000000000000
            2823	0.0500000000000000	0	0	0.0500000000000000
            2824	0.0500000000000000	0	0	0.0500000000000000
            2825	0.0500000000000000	0	0	0.0500000000000000
            2826	0.0500000000000000	0	0	0.0500000000000000
            2827	0.0500000000000000	0	0	0.0500000000000000
            2828	0.0500000000000000	0	0	0.0500000000000000
            2829	0.0500000000000000	0	0	0.0500000000000000
            2830	0.0500000000000000	0	0	0.0500000000000000
            2831	0.0500000000000000	0	0	0.0500000000000000
            2832	0.0500000000000000	0	0	0.0500000000000000
            2833	0.0500000000000000	0	0	0.0500000000000000
            2834	0.0500000000000000	0	0	0.0500000000000000
            2835	0.0500000000000000	0	0	0.0500000000000000
            2836	0.0500000000000000	0	0	0.0500000000000000
            2837	0.0500000000000000	0	0	0.0500000000000000
            2838	0.0500000000000000	0	0	0.0500000000000000
            2839	0.0500000000000000	0	0	0.0500000000000000
            2840	0.0500000000000000	0	0	0.0500000000000000
            2841	0.0500000000000000	0	0	0.0500000000000000
            2842	0.0500000000000000	0	0	0.0500000000000000
            2843	0.0500000000000000	0	0	0.0500000000000000
            2844	0.0500000000000000	0	0	0.0500000000000000
            2845	0.0500000000000000	0	0	0.0500000000000000
            2846	0.0500000000000000	0	0	0.0500000000000000
            2847	0.0500000000000000	0	0	0.0500000000000000
            2848	0.0500000000000000	0	0	0.0500000000000000
            2849	0.0500000000000000	0	0	0.0500000000000000
            2850	0.0500000000000000	0	0	0.0500000000000000
            2851	0.0500000000000000	0	0	0.0500000000000000
            2852	0.0500000000000000	0	0	0.0500000000000000
            2853	0.0500000000000000	0	0	0.0500000000000000
            2854	0.0500000000000000	0	0	0.0500000000000000
            2855	0.0500000000000000	0	0	0.0500000000000000
            2856	0.0500000000000000	0	0	0.0500000000000000
            2857	0.0500000000000000	0	0	0.0500000000000000
            2858	0.0500000000000000	0	0	0.0500000000000000
            2859	0.0500000000000000	0	0	0.0500000000000000
            2860	0.0500000000000000	0	0	0.0500000000000000
            2861	0.0500000000000000	0	0	0.0500000000000000
            2862	0.0500000000000000	0	0	0.0500000000000000
            2863	0.0500000000000000	0	0	0.0500000000000000
            2864	0.0500000000000000	0	0	0.0500000000000000
            2865	0.0500000000000000	0	0	0.0500000000000000
            2866	0.0500000000000000	0	0	0.0500000000000000
            2867	0.0500000000000000	0	0	0.0500000000000000
            2868	0.0500000000000000	0	0	0.0500000000000000
            2869	0.0933141320853138	0	0	0.0933141320853138
            2870	0.0563419029652256	0	0	0.0563419029652256
            2871	0.123566780948183	0	0	0.123566780948183
            2872	0.190954849702034	0	0	0.190954849702034
            2873	0.326465778364573	0	0	0.326465778364573
            2874	0.230956531711394	0	0	0.230956531711394
            2875	0.371874401433239	0	0	0.371874401433239
            2876	0.470657294164945	0	0	0.470657294164945
            2877	0.583082946912755	0	0	0.583082946912755
            2878	0.525183845245952	0	0	0.525183845245952
            2879	0.672864879599907	0	0	0.672864879599907
            2880	0.646605469965197	0	0	0.646605469965197
            2881	0.693885534604820	0	0	0.693885534604820
            2882	0.826843740017205	0	0	0.826843740017205
            2883	0.950000000000000	0	0	0.950000000000000
            2884	0.791555933538746	0	0	0.791555933538746
            2885	0.950000000000000	0	0	0.950000000000000
            2886	0.950000000000000	0	0	0.950000000000000
            2887	0.950000000000000	0	0	0.950000000000000
            2888	0.950000000000000	0	0	0.950000000000000
            2889	0.950000000000000	0	0	0.950000000000000
            2890	0.950000000000000	0	0	0.950000000000000
            2891	0.950000000000000	0	0	0.950000000000000
            2892	0.950000000000000	0	0	0.950000000000000
            2893	0.950000000000000	0	0	0.950000000000000
            2894	0.950000000000000	0	0	0.950000000000000
            2895	0.950000000000000	0	0	0.950000000000000
            2896	0.760815705678567	0	0	0.760815705678567
            2897	0.435566644513477	0	0	0.435566644513477
            2898	0.570994968722570	0	0	0.570994968722570
            2899	0.270955074852489	0	0	0.270955074852489
            2900	0.207034953098845	0	0	0.207034953098845
            2901	0.0816090615116449	0	0	0.0816090615116449
            2902	0.106656162268603	0	0	0.106656162268603
            2903	0.0500000000000000	0	0	0.0500000000000000
            2904	0.0500000000000000	0	0	0.0500000000000000
            2905	0.0500000000000000	0	0	0.0500000000000000
            2906	0.0500000000000000	0	0	0.0500000000000000
            2907	0.0500000000000000	0	0	0.0500000000000000
            2908	0.0500000000000000	0	0	0.0500000000000000
            2909	0.0500000000000000	0	0	0.0500000000000000
            2910	0.0500000000000000	0	0	0.0500000000000000
            2911	0.0500000000000000	0	0	0.0500000000000000
            2912	0.0500000000000000	0	0	0.0500000000000000
            2913	0.0500000000000000	0	0	0.0500000000000000
            2914	0.0500000000000000	0	0	0.0500000000000000
            2915	0.0500000000000000	0	0	0.0500000000000000
            2916	0.0500000000000000	0	0	0.0500000000000000
            2917	0.0500000000000000	0	0	0.0500000000000000
            2918	0.0500000000000000	0	0	0.0500000000000000
            2919	0.0500000000000000	0	0	0.0500000000000000
            2920	0.0500000000000000	0	0	0.0500000000000000
            2921	0.0500000000000000	0	0	0.0500000000000000
            2922	0.0500000000000000	0	0	0.0500000000000000
            2923	0.0500000000000000	0	0	0.0500000000000000
            2924	0.0500000000000000	0	0	0.0500000000000000
            2925	0.0596201906681783	0	0	0.0596201906681783
            2926	0.0500000000000000	0	0	0.0500000000000000
            2927	0.0932699669065373	0	0	0.0932699669065373
            2928	0.143581168993861	0	0	0.143581168993861
            2929	0.284432014741379	0	0	0.284432014741379
            2930	0.184766318088488	0	0	0.184766318088488
            2931	0.301079132557656	0	0	0.301079132557656
            2932	0.463485690052989	0	0	0.463485690052989
            2933	0.621256637118924	0	0	0.621256637118924
            2934	0.403566740861722	0	0	0.403566740861722
            2935	0.444966113078214	0	0	0.444966113078214
            2936	0.684987445833282	0	0	0.684987445833282
            2937	0.621256637117175	0	0	0.621256637117175
            2938	0.403566740855419	0	0	0.403566740855419
            2939	0.301079132502257	0	0	0.301079132502257
            2940	0.463485690044859	0	0	0.463485690044859
            2941	0.284432014692513	0	0	0.284432014692513
            2942	0.184766317663870	0	0	0.184766317663870
            2943	0.0932699641962764	0	0	0.0932699641962764
            2944	0.143581168700395	0	0	0.143581168700395
            2945	0.0500000000000000	0	0	0.0500000000000000
            2946	0.0500000000000000	0	0	0.0500000000000000
            2947	0.0500000000000000	0	0	0.0500000000000000
            2948	0.0500000000000000	0	0	0.0500000000000000
            2949	0.0500000000000000	0	0	0.0500000000000000
            2950	0.0500000000000000	0	0	0.0500000000000000
            2951	0.0500000000000000	0	0	0.0500000000000000
            2952	0.0500000000000000	0	0	0.0500000000000000
            2953	0.0500000000000000	0	0	0.0500000000000000
            2954	0.0500000000000000	0	0	0.0500000000000000
            2955	0.0500000000000000	0	0	0.0500000000000000
            2956	0.0500000000000000	0	0	0.0500000000000000
            2957	0.0500000000000000	0	0	0.0500000000000000
            2958	0.0500000000000000	0	0	0.0500000000000000
            2959	0.0500000000000000	0	0	0.0500000000000000
            2960	0.0500000000000000	0	0	0.0500000000000000
            2961	0.0500000000000000	0	0	0.0500000000000000
            2962	0.0500000000000000	0	0	0.0500000000000000
            2963	0.0500000000000000	0	0	0.0500000000000000
            2964	0.0500000000000000	0	0	0.0500000000000000
            2965	0.0500000000000000	0	0	0.0500000000000000
            2966	0.0500000000000000	0	0	0.0500000000000000
            2967	0.0500000000000000	0	0	0.0500000000000000
            2968	0.0500000000000000	0	0	0.0500000000000000
            2969	0.0500000000000000	0	0	0.0500000000000000
            2970	0.0500000000000000	0	0	0.0500000000000000
            2971	0.0500000000000000	0	0	0.0500000000000000
            2972	0.0500000000000000	0	0	0.0500000000000000
            2973	0.0500000000000000	0	0	0.0500000000000000
            2974	0.0500000000000000	0	0	0.0500000000000000
            2975	0.0500000000000000	0	0	0.0500000000000000
            2976	0.0500000000000000	0	0	0.0500000000000000
            2977	0.0500000000000000	0	0	0.0500000000000000
            2978	0.0500000000000000	0	0	0.0500000000000000
            2979	0.0500000000000000	0	0	0.0500000000000000
            2980	0.0500000000000000	0	0	0.0500000000000000
            2981	0.0500000000000000	0	0	0.0500000000000000
            2982	0.0500000000000000	0	0	0.0500000000000000
            2983	0.0500000000000000	0	0	0.0500000000000000
            2984	0.0500000000000000	0	0	0.0500000000000000
            2985	0.0500000000000000	0	0	0.0500000000000000
            2986	0.0500000000000000	0	0	0.0500000000000000
            2987	0.0500000000000000	0	0	0.0500000000000000
            2988	0.0500000000000000	0	0	0.0500000000000000
            2989	0.0500000000000000	0	0	0.0500000000000000
            2990	0.0500000000000000	0	0	0.0500000000000000
            2991	0.0500000000000000	0	0	0.0500000000000000
            2992	0.0500000000000000	0	0	0.0500000000000000
            2993	0.0500000000000000	0	0	0.0500000000000000
            2994	0.0500000000000000	0	0	0.0500000000000000
            2995	0.0500000000000000	0	0	0.0500000000000000
            2996	0.0500000000000000	0	0	0.0500000000000000
            2997	0.0500000000000000	0	0	0.0500000000000000
            2998	0.0500000000000000	0	0	0.0500000000000000
            2999	0.0917802336447638	0	0	0.0917802336447638
            3000	0.0954732898829486	0	0	0.0954732898829486
            3001	0.203033178924346	0	0	0.203033178924346
            3002	0.212400520718594	0	0	0.212400520718594
            3003	0.412957117223757	0	0	0.412957117223757
            3004	0.370198771177609	0	0	0.370198771177609
            3005	0.582974530628741	0	0	0.582974530628741
            3006	0.677987629830281	0	0	0.677987629830281
            3007	0.950000000000000	0	0	0.950000000000000
            3008	0.809447802236943	0	0	0.809447802236943
            3009	0.950000000000000	0	0	0.950000000000000
            3010	0.950000000000000	0	0	0.950000000000000
            3011	0.950000000000000	0	0	0.950000000000000
            3012	0.950000000000000	0	0	0.950000000000000
            3013	0.950000000000000	0	0	0.950000000000000
            3014	0.950000000000000	0	0	0.950000000000000
            3015	0.950000000000000	0	0	0.950000000000000
            3016	0.950000000000000	0	0	0.950000000000000
            3017	0.950000000000000	0	0	0.950000000000000
            3018	0.950000000000000	0	0	0.950000000000000
            3019	0.950000000000000	0	0	0.950000000000000
            3020	0.950000000000000	0	0	0.950000000000000
            3021	0.950000000000000	0	0	0.950000000000000
            3022	0.950000000000000	0	0	0.950000000000000
            3023	0.950000000000000	0	0	0.950000000000000
            3024	0.950000000000000	0	0	0.950000000000000
            3025	0.626962723608423	0	0	0.626962723608423
            3026	0.568628577297933	0	0	0.568628577297933
            3027	0.269491586076553	0	0	0.269491586076553
            3028	0.297308360684687	0	0	0.297308360684687
            3029	0.116986538270913	0	0	0.116986538270913
            3030	0.106064544323192	0	0	0.106064544323192
            3031	0.0500000000000000	0	0	0.0500000000000000
            3032	0.0500000000000000	0	0	0.0500000000000000
            3033	0.0500000000000000	0	0	0.0500000000000000
            3034	0.0500000000000000	0	0	0.0500000000000000
            3035	0.0500000000000000	0	0	0.0500000000000000
            3036	0.0500000000000000	0	0	0.0500000000000000
            3037	0.0500000000000000	0	0	0.0500000000000000
            3038	0.0500000000000000	0	0	0.0500000000000000
            3039	0.0500000000000000	0	0	0.0500000000000000
            3040	0.0500000000000000	0	0	0.0500000000000000
            3041	0.0500000000000000	0	0	0.0500000000000000
            3042	0.0500000000000000	0	0	0.0500000000000000
            3043	0.0500000000000000	0	0	0.0500000000000000
            3044	0.0500000000000000	0	0	0.0500000000000000
            3045	0.0500000000000000	0	0	0.0500000000000000
            3046	0.0500000000000000	0	0	0.0500000000000000
            3047	0.0500000000000000	0	0	0.0500000000000000
            3048	0.0500000000000000	0	0	0.0500000000000000
            3049	0.0500000000000000	0	0	0.0500000000000000
            3050	0.0500000000000000	0	0	0.0500000000000000
            3051	0.0500000000000000	0	0	0.0500000000000000
            3052	0.0500000000000000	0	0	0.0500000000000000
            3053	0.0500000000000000	0	0	0.0500000000000000
            3054	0.0500000000000000	0	0	0.0500000000000000
            3055	0.0621813023718655	0	0	0.0621813023718655
            3056	0.0621811887469732	0	0	0.0621811887469732
            3057	0.123179908851830	0	0	0.123179908851830
            3058	0.123179926621328	0	0	0.123179926621328
            3059	0.200723270762132	0	0	0.200723270762132
            3060	0.200723268473953	0	0	0.200723268473953
            3061	0.269049649642224	0	0	0.269049649642224
            3062	0.269049649889431	0	0	0.269049649889431
            3063	0.296649760002909	0	0	0.296649760002909
            3064	0.296649759972835	0	0	0.296649759972835
            3065	0.269049649597753	0	0	0.269049649597753
            3066	0.269049649611536	0	0	0.269049649611536
            3067	0.200723268064124	0	0	0.200723268064124
            3068	0.200723268050426	0	0	0.200723268050426
            3069	0.123179905558281	0	0	0.123179905558281
            3070	0.123179905570635	0	0	0.123179905570635
            3071	0.0621811676939353	0	0	0.0621811676939353
            3072	0.0621811676847427	0	0	0.0621811676847427
            3073	0.0500000000000000	0	0	0.0500000000000000
            3074	0.0500000000000000	0	0	0.0500000000000000
            3075	0.0500000000000000	0	0	0.0500000000000000
            3076	0.0500000000000000	0	0	0.0500000000000000
            3077	0.0500000000000000	0	0	0.0500000000000000
            3078	0.0500000000000000	0	0	0.0500000000000000
            3079	0.0500000000000000	0	0	0.0500000000000000
            3080	0.0500000000000000	0	0	0.0500000000000000
            3081	0.0500000000000000	0	0	0.0500000000000000
            3082	0.0500000000000000	0	0	0.0500000000000000
            3083	0.0500000000000000	0	0	0.0500000000000000
            3084	0.0500000000000000	0	0	0.0500000000000000
            3085	0.0500000000000000	0	0	0.0500000000000000
            3086	0.0500000000000000	0	0	0.0500000000000000
            3087	0.0500000000000000	0	0	0.0500000000000000
            3088	0.0500000000000000	0	0	0.0500000000000000
            3089	0.0500000000000000	0	0	0.0500000000000000
            3090	0.0500000000000000	0	0	0.0500000000000000
            3091	0.0500000000000000	0	0	0.0500000000000000
            3092	0.0500000000000000	0	0	0.0500000000000000
            3093	0.0500000000000000	0	0	0.0500000000000000
            3094	0.0500000000000000	0	0	0.0500000000000000
            3095	0.0500000000000000	0	0	0.0500000000000000
            3096	0.0500000000000000	0	0	0.0500000000000000
            3097	0.0500000000000000	0	0	0.0500000000000000
            3098	0.0500000000000000	0	0	0.0500000000000000
            3099	0.0500000000000000	0	0	0.0500000000000000
            3100	0.0500000000000000	0	0	0.0500000000000000
            3101	0.0500000000000000	0	0	0.0500000000000000
            3102	0.0500000000000000	0	0	0.0500000000000000
            3103	0.0500000000000000	0	0	0.0500000000000000
            3104	0.0500000000000000	0	0	0.0500000000000000
            3105	0.0500000000000000	0	0	0.0500000000000000
            3106	0.0500000000000000	0	0	0.0500000000000000
            3107	0.0500000000000000	0	0	0.0500000000000000
            3108	0.0500000000000000	0	0	0.0500000000000000
            3109	0.0500000000000000	0	0	0.0500000000000000
            3110	0.0500000000000000	0	0	0.0500000000000000
            3111	0.0500000000000000	0	0	0.0500000000000000
            3112	0.0500000000000000	0	0	0.0500000000000000
            3113	0.0500000000000000	0	0	0.0500000000000000
            3114	0.0500000000000000	0	0	0.0500000000000000
            3115	0.0500000000000000	0	0	0.0500000000000000
            3116	0.0500000000000000	0	0	0.0500000000000000
            3117	0.0500000000000000	0	0	0.0500000000000000
            3118	0.0500000000000000	0	0	0.0500000000000000
            3119	0.0500000000000000	0	0	0.0500000000000000
            3120	0.0500000000000000	0	0	0.0500000000000000
            3121	0.0500000000000000	0	0	0.0500000000000000
            3122	0.0500000000000000	0	0	0.0500000000000000
            3123	0.0500000000000000	0	0	0.0500000000000000
            3124	0.0500000000000000	0	0	0.0500000000000000
            3125	0.0500000000000000	0	0	0.0500000000000000
            3126	0.0500000000000000	0	0	0.0500000000000000
            3127	0.0836619371415612	0	0	0.0836619371415612
            3128	0.0915016987038106	0	0	0.0915016987038106
            3129	0.217405937529498	0	0	0.217405937529498
            3130	0.195481652689385	0	0	0.195481652689385
            3131	0.384772491254089	0	0	0.384772491254089
            3132	0.430176499839057	0	0	0.430176499839057
            3133	0.710987188829704	0	0	0.710987188829704
            3134	0.632654544632007	0	0	0.632654544632007
            3135	0.873399747048543	0	0	0.873399747048543
            3136	0.950000000000000	0	0	0.950000000000000
            3137	0.950000000000000	0	0	0.950000000000000
            3138	0.950000000000000	0	0	0.950000000000000
            3139	0.950000000000000	0	0	0.950000000000000
            3140	0.950000000000000	0	0	0.950000000000000
            3141	0.950000000000000	0	0	0.950000000000000
            3142	0.950000000000000	0	0	0.950000000000000
            3143	0.950000000000000	0	0	0.950000000000000
            3144	0.950000000000000	0	0	0.950000000000000
            3145	0.950000000000000	0	0	0.950000000000000
            3146	0.950000000000000	0	0	0.950000000000000
            3147	0.870807352083124	0	0	0.870807352083124
            3148	0.950000000000000	0	0	0.950000000000000
            3149	0.950000000000000	0	0	0.950000000000000
            3150	0.684556534789910	0	0	0.684556534789910
            3151	0.461036068742155	0	0	0.461036068742155
            3152	0.747207495880899	0	0	0.747207495880899
            3153	0.424805606786689	0	0	0.424805606786689
            3154	0.261462336573664	0	0	0.261462336573664
            3155	0.124241332847986	0	0	0.124241332847986
            3156	0.201270376818960	0	0	0.201270376818960
            3157	0.0794759372909423	0	0	0.0794759372909423
            3158	0.0502897425077268	0	0	0.0502897425077268
            3159	0.0500000000000000	0	0	0.0500000000000000
            3160	0.0500000000000000	0	0	0.0500000000000000
            3161	0.0500000000000000	0	0	0.0500000000000000
            3162	0.0500000000000000	0	0	0.0500000000000000
            3163	0.0500000000000000	0	0	0.0500000000000000
            3164	0.0500000000000000	0	0	0.0500000000000000
            3165	0.0500000000000000	0	0	0.0500000000000000
            3166	0.0500000000000000	0	0	0.0500000000000000
            3167	0.0500000000000000	0	0	0.0500000000000000
            3168	0.0500000000000000	0	0	0.0500000000000000
            3169	0.0500000000000000	0	0	0.0500000000000000
            3170	0.0500000000000000	0	0	0.0500000000000000
            3171	0.0500000000000000	0	0	0.0500000000000000
            3172	0.0500000000000000	0	0	0.0500000000000000
            3173	0.0500000000000000	0	0	0.0500000000000000
            3174	0.0500000000000000	0	0	0.0500000000000000
            3175	0.0500000000000000	0	0	0.0500000000000000
            3176	0.0500000000000000	0	0	0.0500000000000000
            3177	0.0500000000000000	0	0	0.0500000000000000
            3178	0.0500000000000000	0	0	0.0500000000000000
            3179	0.0500000000000000	0	0	0.0500000000000000
            3180	0.0500000000000000	0	0	0.0500000000000000
            3181	0.0500000000000000	0	0	0.0500000000000000
            3182	0.0596363139235384	0	0	0.0596363139235384
            3183	0.143584234306817	0	0	0.143584234306817
            3184	0.0932706726763710	0	0	0.0932706726763710
            3185	0.184766428479483	0	0	0.184766428479483
            3186	0.284432494492251	0	0	0.284432494492251
            3187	0.463485752729726	0	0	0.463485752729726
            3188	0.301079146813963	0	0	0.301079146813963
            3189	0.403566742482475	0	0	0.403566742482475
            3190	0.621256645663420	0	0	0.621256645663420
            3191	0.684987449694736	0	0	0.684987449694736
            3192	0.444966113402925	0	0	0.444966113402925
            3193	0.403566741125256	0	0	0.403566741125256
            3194	0.621256641573298	0	0	0.621256641573298
            3195	0.463485694915965	0	0	0.463485694915965
            3196	0.301079132790209	0	0	0.301079132790209
            3197	0.184766317924585	0	0	0.184766317924585
            3198	0.284432019108406	0	0	0.284432019108406
            3199	0.143581171994709	0	0	0.143581171994709
            3200	0.0932699643907109	0	0	0.0932699643907109
            3201	0.0500000000000000	0	0	0.0500000000000000
            3202	0.0500000000000000	0	0	0.0500000000000000
            3203	0.0500000000000000	0	0	0.0500000000000000
            3204	0.0500000000000000	0	0	0.0500000000000000
            3205	0.0500000000000000	0	0	0.0500000000000000
            3206	0.0500000000000000	0	0	0.0500000000000000
            3207	0.0500000000000000	0	0	0.0500000000000000
            3208	0.0500000000000000	0	0	0.0500000000000000
            3209	0.0500000000000000	0	0	0.0500000000000000
            3210	0.0500000000000000	0	0	0.0500000000000000
            3211	0.0500000000000000	0	0	0.0500000000000000
            3212	0.0500000000000000	0	0	0.0500000000000000
            3213	0.0500000000000000	0	0	0.0500000000000000
            3214	0.0500000000000000	0	0	0.0500000000000000
            3215	0.0500000000000000	0	0	0.0500000000000000
            3216	0.0500000000000000	0	0	0.0500000000000000
            3217	0.0500000000000000	0	0	0.0500000000000000
            3218	0.0500000000000000	0	0	0.0500000000000000
            3219	0.0500000000000000	0	0	0.0500000000000000
            3220	0.0500000000000000	0	0	0.0500000000000000
            3221	0.0500000000000000	0	0	0.0500000000000000
            3222	0.0614806887411965	0	0	0.0614806887411965
            3223	0.0865160496095424	0	0	0.0865160496095424
            3224	0.0500000000000000	0	0	0.0500000000000000
            3225	0.0500000000000000	0	0	0.0500000000000000
            3226	0.105282601205030	0	0	0.105282601205030
            3227	0.115081910287069	0	0	0.115081910287069
            3228	0.0500000000000000	0	0	0.0500000000000000
            3229	0.0500000000000000	0	0	0.0500000000000000
            3230	0.117786239323911	0	0	0.117786239323911
            3231	0.115082172952402	0	0	0.115082172952402
            3232	0.0500000000000000	0	0	0.0500000000000000
            3233	0.0500000000000000	0	0	0.0500000000000000
            3234	0.105284675374578	0	0	0.105284675374578
            3235	0.0865294853344717	0	0	0.0865294853344717
            3236	0.0500000000000000	0	0	0.0500000000000000
            3237	0.0500000000000000	0	0	0.0500000000000000
            3238	0.0615527653950519	0	0	0.0615527653950519
            3239	0.0500000000000000	0	0	0.0500000000000000
            3240	0.0500000000000000	0	0	0.0500000000000000
            3241	0.0500000000000000	0	0	0.0500000000000000
            3242	0.0500000000000000	0	0	0.0500000000000000
            3243	0.0500000000000000	0	0	0.0500000000000000
            3244	0.0500000000000000	0	0	0.0500000000000000
            3245	0.0500000000000000	0	0	0.0500000000000000
            3246	0.0500000000000000	0	0	0.0500000000000000
            3247	0.0500000000000000	0	0	0.0500000000000000
            3248	0.0500000000000000	0	0	0.0500000000000000
            3249	0.0500000000000000	0	0	0.0500000000000000
            3250	0.0500000000000000	0	0	0.0500000000000000
            3251	0.0522582827144322	0	0	0.0522582827144322
            3252	0.0500000000000000	0	0	0.0500000000000000
            3253	0.0500000000000000	0	0	0.0500000000000000
            3254	0.0711434714759795	0	0	0.0711434714759795
            3255	0.0950741689381561	0	0	0.0950741689381561
            3256	0.0754689981614220	0	0	0.0754689981614220
            3257	0.156144237327163	0	0	0.156144237327163
            3258	0.136259493839931	0	0	0.136259493839931
            3259	0.207212317448973	0	0	0.207212317448973
            3260	0.292832189196074	0	0	0.292832189196074
            3261	0.471559423560837	0	0	0.471559423560837
            3262	0.303745068806427	0	0	0.303745068806427
            3263	0.396053484430067	0	0	0.396053484430067
            3264	0.641045436291754	0	0	0.641045436291754
            3265	0.740071164304380	0	0	0.740071164304380
            3266	0.444738075878975	0	0	0.444738075878975
            3267	0.432203318848242	0	0	0.432203318848242
            3268	0.745137222663026	0	0	0.745137222663026
            3269	0.686868143735415	0	0	0.686868143735415
            3270	0.376332032755298	0	0	0.376332032755298
            3271	0.310332833990324	0	0	0.310332833990324
            3272	0.612955317581832	0	0	0.612955317581832
            3273	0.541567824594360	0	0	0.541567824594360
            3274	0.253037197991306	0	0	0.253037197991306
            3275	0.202154802070828	0	0	0.202154802070828
            3276	0.458629085130277	0	0	0.458629085130277
            3277	0.352114105058240	0	0	0.352114105058240
            3278	0.150346595719145	0	0	0.150346595719145
            3279	0.0992301486135908	0	0	0.0992301486135908
            3280	0.234769319834151	0	0	0.234769319834151
            3281	0.133165684252971	0	0	0.133165684252971
            3282	0.0577749806825027	0	0	0.0577749806825027
            3283	0.0500000000000000	0	0	0.0500000000000000
            3284	0.0650172383111223	0	0	0.0650172383111223
            3285	0.0500000000000000	0	0	0.0500000000000000
            3286	0.0500000000000000	0	0	0.0500000000000000
            3287	0.0500000000000000	0	0	0.0500000000000000
            3288	0.0500000000000000	0	0	0.0500000000000000
            3289	0.0500000000000000	0	0	0.0500000000000000
            3290	0.0597884699547134	0	0	0.0597884699547134
            3291	0.0791403713607206	0	0	0.0791403713607206
            3292	0.0500000000000000	0	0	0.0500000000000000
            3293	0.0500000000000000	0	0	0.0500000000000000
            3294	0.0870753757720253	0	0	0.0870753757720253
            3295	0.0789467910788849	0	0	0.0789467910788849
            3296	0.0500000000000000	0	0	0.0500000000000000
            3297	0.0500000000000000	0	0	0.0500000000000000
            3298	0.0589012015973186	0	0	0.0589012015973186
            3299	0.0500000000000000	0	0	0.0500000000000000
            3300	0.0500000000000000	0	0	0.0500000000000000
            3301	0.0500000000000000	0	0	0.0500000000000000
            3302	0.0500000000000000	0	0	0.0500000000000000
            3303	0.0500000000000000	0	0	0.0500000000000000
            3304	0.0500000000000000	0	0	0.0500000000000000
            3305	0.0500000000000000	0	0	0.0500000000000000
            3306	0.0500000000000000	0	0	0.0500000000000000
            3307	0.0500000000000000	0	0	0.0500000000000000
            3308	0.0500000000000000	0	0	0.0500000000000000
            3309	0.0791664292978050	0	0	0.0791664292978050
            3310	0.0872433358063170	0	0	0.0872433358063170
            3311	0.209732309795269	0	0	0.209732309795269
            3312	0.190526378379810	0	0	0.190526378379810
            3313	0.377409823927947	0	0	0.377409823927947
            3314	0.415418113456686	0	0	0.415418113456686
            3315	0.676921807820469	0	0	0.676921807820469
            3316	0.614991690739599	0	0	0.614991690739599
            3317	0.824335155313315	0	0	0.824335155313315
            3318	0.907345913551570	0	0	0.907345913551570
            3319	0.950000000000000	0	0	0.950000000000000
            3320	0.908898502491208	0	0	0.908898502491208
            3321	0.824335165667583	0	0	0.824335165667583
            3322	0.907346230094840	0	0	0.907346230094840
            3323	0.676921782701030	0	0	0.676921782701030
            3324	0.614991525012735	0	0	0.614991525012735
            3325	0.377408173488937	0	0	0.377408173488937
            3326	0.415413758729962	0	0	0.415413758729962
            3327	0.209700872283885	0	0	0.209700872283885
            3328	0.190515511822688	0	0	0.190515511822688
            3329	0.0500000000000000	0	0	0.0500000000000000
            3330	0.0500000000000000	0	0	0.0500000000000000
            3331	0.0500000000000000	0	0	0.0500000000000000
            3332	0.0500000000000000	0	0	0.0500000000000000
            3333	0.0500000000000000	0	0	0.0500000000000000
            3334	0.0500000000000000	0	0	0.0500000000000000
            3335	0.0500000000000000	0	0	0.0500000000000000
            3336	0.0500000000000000	0	0	0.0500000000000000
            3337	0.0500000000000000	0	0	0.0500000000000000
            3338	0.0500000000000000	0	0	0.0500000000000000
            3339	0.0500000000000000	0	0	0.0500000000000000
            3340	0.0500000000000000	0	0	0.0500000000000000
            3341	0.0500000000000000	0	0	0.0500000000000000
            3342	0.0500000000000000	0	0	0.0500000000000000
            3343	0.0500000000000000	0	0	0.0500000000000000
            3344	0.0500000000000000	0	0	0.0500000000000000
            3345	0.0500000000000000	0	0	0.0500000000000000
            3346	0.0878406717820479	0	0	0.0878406717820479
            3347	0.175890156454177	0	0	0.175890156454177
            3348	0.0887892346117761	0	0	0.0887892346117761
            3349	0.148061744686899	0	0	0.148061744686899
            3350	0.293307979059191	0	0	0.293307979059191
            3351	0.412744989284262	0	0	0.412744989284262
            3352	0.208353506664610	0	0	0.208353506664610
            3353	0.253548321419834	0	0	0.253548321419834
            3354	0.502275198282168	0	0	0.502275198282168
            3355	0.549025096571479	0	0	0.549025096571479
            3356	0.277147648773793	0	0	0.277147648773793
            3357	0.283660406916031	0	0	0.283660406916031
            3358	0.561926826594981	0	0	0.561926826594981
            3359	0.549027288385026	0	0	0.549027288385026
            3360	0.277148468733995	0	0	0.277148468733995
            3361	0.253554728683791	0	0	0.253554728683791
            3362	0.502292058502008	0	0	0.502292058502008
            3363	0.412851773996746	0	0	0.412851773996746
            3364	0.208394650328225	0	0	0.208394650328225
            3365	0.148280899011444	0	0	0.148280899011444
            3366	0.293870210500093	0	0	0.293870210500093
            3367	0.178351232372676	0	0	0.178351232372676
            3368	0.0897577827265852	0	0	0.0897577827265852
            3369	0.0500000000000000	0	0	0.0500000000000000
            3370	0.0967977979449579	0	0	0.0967977979449579
            3371	0.0634050568202221	0	0	0.0634050568202221
            3372	0.0500000000000000	0	0	0.0500000000000000
            3373	0.0500000000000000	0	0	0.0500000000000000
            3374	0.0806106400444504	0	0	0.0806106400444504
            3375	0.146362096593463	0	0	0.146362096593463
            3376	0.0592722091719676	0	0	0.0592722091719676
            3377	0.101063583773783	0	0	0.101063583773783
            3378	0.249758843628260	0	0	0.249758843628260
            3379	0.361287615706735	0	0	0.361287615706735
            3380	0.147244403216645	0	0	0.147244403216645
            3381	0.181552348880002	0	0	0.181552348880002
            3382	0.436871686340243	0	0	0.436871686340243
            3383	0.442684288686459	0	0	0.442684288686459
            3384	0.193899231161876	0	0	0.193899231161876
            3385	0.190221144281526	0	0	0.190221144281526
            3386	0.380635334721229	0	0	0.380635334721229
            3387	0.286963333419118	0	0	0.286963333419118
            3388	0.188441877257521	0	0	0.188441877257521
            3389	0.200991881836072	0	0	0.200991881836072
            3390	0.203258963475332	0	0	0.203258963475332
            3391	0.148577490503868	0	0	0.148577490503868
            3392	0.221215525017320	0	0	0.221215525017320
            3393	0.229341948892721	0	0	0.229341948892721
            3394	0.116649927134743	0	0	0.116649927134743
            3395	0.0930423473145340	0	0	0.0930423473145340
            3396	0.212255249316287	0	0	0.212255249316287
            3397	0.175257849164580	0	0	0.175257849164580
            3398	0.0704206705386706	0	0	0.0704206705386706
            3399	0.0500000000000000	0	0	0.0500000000000000
            3400	0.134335785381105	0	0	0.134335785381105
            3401	0.100681649909101	0	0	0.100681649909101
            3402	0.0500000000000000	0	0	0.0500000000000000
            3403	0.0500000000000000	0	0	0.0500000000000000
            3404	0.0750197898098735	0	0	0.0750197898098735
            3405	0.0537160409531194	0	0	0.0537160409531194
            3406	0.0500000000000000	0	0	0.0500000000000000
            3407	0.0500000000000000	0	0	0.0500000000000000
            3408	0.0500000000000000	0	0	0.0500000000000000
            3409	0.0500000000000000	0	0	0.0500000000000000
            3410	0.0500000000000000	0	0	0.0500000000000000
            3411	0.0500000000000000	0	0	0.0500000000000000
            3412	0.0500000000000000	0	0	0.0500000000000000
            3413	0.0500000000000000	0	0	0.0500000000000000
            3414	0.0880187014016519	0	0	0.0880187014016519
            3415	0.172740358797985	0	0	0.172740358797985
            3416	0.0881720643160812	0	0	0.0881720643160812
            3417	0.142136036765325	0	0	0.142136036765325
            3418	0.281049254957726	0	0	0.281049254957726
            3419	0.376622650128651	0	0	0.376622650128651
            3420	0.190177593315466	0	0	0.190177593315466
            3421	0.209624066594621	0	0	0.209624066594621
            3422	0.415240420080678	0	0	0.415240420080678
            3423	0.376604521843081	0	0	0.376604521843081
            3424	0.190111907732056	0	0	0.190111907732056
            3425	0.141837100505042	0	0	0.141837100505042
            3426	0.280968370169125	0	0	0.280968370169125
            3427	0.172460574863486	0	0	0.172460574863486
            3428	0.0870904685387948	0	0	0.0870904685387948
            3429	0.0500000000000000	0	0	0.0500000000000000
            3430	0.0872874026182158	0	0	0.0872874026182158
            3431	0.0500000000000000	0	0	0.0500000000000000
            3432	0.0500000000000000	0	0	0.0500000000000000
            3433	0.0500000000000000	0	0	0.0500000000000000
            3434	0.0500000000000000	0	0	0.0500000000000000
            3435	0.0500000000000000	0	0	0.0500000000000000
            3436	0.0500000000000000	0	0	0.0500000000000000
            3437	0.0793820165915252	0	0	0.0793820165915252
            3438	0.0598638466320753	0	0	0.0598638466320753
            3439	0.142390448653131	0	0	0.142390448653131
            3440	0.190275376061611	0	0	0.190275376061611
            3441	0.376792344626665	0	0	0.376792344626665
            3442	0.281796779593685	0	0	0.281796779593685
            3443	0.459159662834297	0	0	0.459159662834297
            3444	0.613970501452071	0	0	0.613970501452071
            3445	0.822965694911208	0	0	0.822965694911208
            3446	0.615461947832715	0	0	0.615461947832715
            3447	0.678610543591698	0	0	0.678610543591698
            3448	0.907390015470775	0	0	0.907390015470775
            3449	0.822969162315434	0	0	0.822969162315434
            3450	0.615489845957342	0	0	0.615489845957342
            3451	0.459201958596789	0	0	0.459201958596789
            3452	0.613974789065385	0	0	0.613974789065385
            3453	0.376786367872965	0	0	0.376786367872965
            3454	0.281820004855245	0	0	0.281820004855245
            3455	0.142275184317680	0	0	0.142275184317680
            3456	0.190203224997470	0	0	0.190203224997470
            3457	0.0500000000000000	0	0	0.0500000000000000
            3458	0.0500000000000000	0	0	0.0500000000000000
            3459	0.0500000000000000	0	0	0.0500000000000000
            3460	0.0500000000000000	0	0	0.0500000000000000
            3461	0.0500000000000000	0	0	0.0500000000000000
            3462	0.0500000000000000	0	0	0.0500000000000000
            3463	0.0500000000000000	0	0	0.0500000000000000
            3464	0.0500000000000000	0	0	0.0500000000000000
            3465	0.0500000000000000	0	0	0.0500000000000000
            3466	0.0500000000000000	0	0	0.0500000000000000
            3467	0.0500000000000000	0	0	0.0500000000000000
            3468	0.0500000000000000	0	0	0.0500000000000000
            3469	0.0500000000000000	0	0	0.0500000000000000
            3470	0.0500000000000000	0	0	0.0500000000000000
            3471	0.0793429933243963	0	0	0.0793429933243963
            3472	0.0591575529536571	0	0	0.0591575529536571
            3473	0.143140096348517	0	0	0.143140096348517
            3474	0.191876243817418	0	0	0.191876243817418
            3475	0.384182193619629	0	0	0.384182193619629
            3476	0.286615800652284	0	0	0.286615800652284
            3477	0.477949266604352	0	0	0.477949266604352
            3478	0.640643959073242	0	0	0.640643959073242
            3479	0.901517996403331	0	0	0.901517996403331
            3480	0.672573371191990	0	0	0.672573371191990
            3481	0.818464010944011	0	0	0.818464010944011
            3482	0.950000000000000	0	0	0.950000000000000
            3483	0.950000000000000	0	0	0.950000000000000
            3484	0.894643597689017	0	0	0.894643597689017
            3485	0.915667301538965	0	0	0.915667301538965
            3486	0.950000000000000	0	0	0.950000000000000
            3487	0.950000000000000	0	0	0.950000000000000
            3488	0.894648707657507	0	0	0.894648707657507
            3489	0.818502454617099	0	0	0.818502454617099
            3490	0.950000000000000	0	0	0.950000000000000
            3491	0.901985142838483	0	0	0.901985142838483
            3492	0.672811964353227	0	0	0.672811964353227
            3493	0.479183407952865	0	0	0.479183407952865
            3494	0.643001955547034	0	0	0.643001955547034
            3495	0.394129366686747	0	0	0.394129366686747
            3496	0.291936937610597	0	0	0.291936937610597
            3497	0.162264054480523	0	0	0.162264054480523
            3498	0.226954242676989	0	0	0.226954242676989
            3499	0.182758000177526	0	0	0.182758000177526
            3500	0.116447013763122	0	0	0.116447013763122
            3501	0.163256579364716	0	0	0.163256579364716
            3502	0.282124792425613	0	0	0.282124792425613
            3503	0.533208929199479	0	0	0.533208929199479
            3504	0.303487724736815	0	0	0.303487724736815
            3505	0.517853814284801	0	0	0.517853814284801
            3506	0.906543073135190	0	0	0.906543073135190
            3507	0.950000000000000	0	0	0.950000000000000
            3508	0.746685778207196	0	0	0.746685778207196
            3509	0.898666681075464	0	0	0.898666681075464
            3510	0.950000000000000	0	0	0.950000000000000
            3511	0.950000000000000	0	0	0.950000000000000
            3512	0.901587866435750	0	0	0.901587866435750
            3513	0.755088372383148	0	0	0.755088372383148
            3514	0.950000000000000	0	0	0.950000000000000
            3515	0.902343284179088	0	0	0.902343284179088
            3516	0.530874228693756	0	0	0.530874228693756
            3517	0.318378490715668	0	0	0.318378490715668
            3518	0.523073281662204	0	0	0.523073281662204
            3519	0.254866316406642	0	0	0.254866316406642
            3520	0.169496667682973	0	0	0.169496667682973
            3521	0.0866121359252040	0	0	0.0866121359252040
            3522	0.106441844045029	0	0	0.106441844045029
            3523	0.0500000000000000	0	0	0.0500000000000000
            3524	0.0500000000000000	0	0	0.0500000000000000
            3525	0.0500000000000000	0	0	0.0500000000000000
            3526	0.0500000000000000	0	0	0.0500000000000000
            3527	0.0500000000000000	0	0	0.0500000000000000
            3528	0.0500000000000000	0	0	0.0500000000000000
            3529	0.0500000000000000	0	0	0.0500000000000000
            3530	0.0500000000000000	0	0	0.0500000000000000
            3531	0.0500000000000000	0	0	0.0500000000000000
            3532	0.0500000000000000	0	0	0.0500000000000000
            3533	0.0500000000000000	0	0	0.0500000000000000
            3534	0.0500000000000000	0	0	0.0500000000000000
            3535	0.0500000000000000	0	0	0.0500000000000000
            3536	0.0500000000000000	0	0	0.0500000000000000
            3537	0.0500000000000000	0	0	0.0500000000000000
            3538	0.0500000000000000	0	0	0.0500000000000000
            3539	0.0790523113187919	0	0	0.0790523113187919
            3540	0.0594766873249532	0	0	0.0594766873249532
            3541	0.142057279890327	0	0	0.142057279890327
            3542	0.190152606225791	0	0	0.190152606225791
            3543	0.376617488391932	0	0	0.376617488391932
            3544	0.281037027568733	0	0	0.281037027568733
            3545	0.457853236503376	0	0	0.457853236503376
            3546	0.613684031724580	0	0	0.613684031724580
            3547	0.822578407290011	0	0	0.822578407290011
            3548	0.613684690407298	0	0	0.613684690407298
            3549	0.676634694369373	0	0	0.676634694369373
            3550	0.906960789950671	0	0	0.906960789950671
            3551	0.822577773142415	0	0	0.822577773142415
            3552	0.613680707138421	0	0	0.613680707138421
            3553	0.457836504687559	0	0	0.457836504687559
            3554	0.613682021243129	0	0	0.613682021243129
            3555	0.376617188503166	0	0	0.376617188503166
            3556	0.280987604615538	0	0	0.280987604615538
            3557	0.141985618451804	0	0	0.141985618451804
            3558	0.190197124848143	0	0	0.190197124848143
            3559	0.0794025186904557	0	0	0.0794025186904557
            3560	0.0597108484371914	0	0	0.0597108484371914
            3561	0.0500000000000000	0	0	0.0500000000000000
            3562	0.0500000000000000	0	0	0.0500000000000000
            3563	0.0500000000000000	0	0	0.0500000000000000
            3564	0.0500000000000000	0	0	0.0500000000000000
            3565	0.0500000000000000	0	0	0.0500000000000000
            3566	0.0500000000000000	0	0	0.0500000000000000
            3567	0.0505587730865215	0	0	0.0505587730865215
            3568	0.0890486105596435	0	0	0.0890486105596435
            3569	0.175961121586398	0	0	0.175961121586398
            3570	0.0996031362731728	0	0	0.0996031362731728
            3571	0.162357339872581	0	0	0.162357339872581
            3572	0.286697542741648	0	0	0.286697542741648
            3573	0.384332755949109	0	0	0.384332755949109
            3574	0.217874867461493	0	0	0.217874867461493
            3575	0.240646943801663	0	0	0.240646943801663
            3576	0.423838797106652	0	0	0.423838797106652
            3577	0.384512691960891	0	0	0.384512691960891
            3578	0.218823413627781	0	0	0.218823413627781
            3579	0.163876453516409	0	0	0.163876453516409
            3580	0.286982618160009	0	0	0.286982618160009
            3581	0.176223198174883	0	0	0.176223198174883
            3582	0.101133907309143	0	0	0.101133907309143
            3583	0.0514747407682362	0	0	0.0514747407682362
            3584	0.0890377000197866	0	0	0.0890377000197866
            3585	0.0500000000000000	0	0	0.0500000000000000
            3586	0.0870390761390196	0	0	0.0870390761390196
            3587	0.0789470675514076	0	0	0.0789470675514076
            3588	0.0500000000000000	0	0	0.0500000000000000
            3589	0.0500000000000000	0	0	0.0500000000000000
            3590	0.0589451611609596	0	0	0.0589451611609596
            3591	0.0500000000000000	0	0	0.0500000000000000
            3592	0.0500000000000000	0	0	0.0500000000000000
            3593	0.0500000000000000	0	0	0.0500000000000000
            3594	0.0500000000000000	0	0	0.0500000000000000
            3595	0.0500000000000000	0	0	0.0500000000000000
            3596	0.0500000000000000	0	0	0.0500000000000000
            3597	0.0500000000000000	0	0	0.0500000000000000
            3598	0.0500000000000000	0	0	0.0500000000000000
            3599	0.0800002394362175	0	0	0.0800002394362175
            3600	0.0876536876390191	0	0	0.0876536876390191
            3601	0.211599232177598	0	0	0.211599232177598
            3602	0.192028142661108	0	0	0.192028142661108
            3603	0.384211071002632	0	0	0.384211071002632
            3604	0.423600557745186	0	0	0.423600557745186
            3605	0.706364768808995	0	0	0.706364768808995
            3606	0.640648475005107	0	0	0.640648475005107
            3607	0.901518578995245	0	0	0.901518578995245
            3608	0.950000000000000	0	0	0.950000000000000
            3609	0.950000000000000	0	0	0.950000000000000
            3610	0.950000000000000	0	0	0.950000000000000
            3611	0.950000000000000	0	0	0.950000000000000
            3612	0.950000000000000	0	0	0.950000000000000
            3613	0.950000000000000	0	0	0.950000000000000
            3614	0.950000000000000	0	0	0.950000000000000
            3615	0.950000000000000	0	0	0.950000000000000
            3616	0.950000000000000	0	0	0.950000000000000
            3617	0.950000000000000	0	0	0.950000000000000
            3618	0.950000000000000	0	0	0.950000000000000
            3619	0.902782750502125	0	0	0.902782750502125
            3620	0.950000000000000	0	0	0.950000000000000
            3621	0.710349997115440	0	0	0.710349997115440
            3622	0.646657651654999	0	0	0.646657651654999
            3623	0.408113276975159	0	0	0.408113276975159
            3624	0.439949348592758	0	0	0.439949348592758
            3625	0.267793379556639	0	0	0.267793379556639
            3626	0.271735489756098	0	0	0.271735489756098
            3627	0.303194972184551	0	0	0.303194972184551
            3628	0.249590578470550	0	0	0.249590578470550
            3629	0.422016266401394	0	0	0.422016266401394
            3630	0.554979784977345	0	0	0.554979784977345
            3631	0.950000000000000	0	0	0.950000000000000
            3632	0.803942160101715	0	0	0.803942160101715
            3633	0.950000000000000	0	0	0.950000000000000
            3634	0.950000000000000	0	0	0.950000000000000
            3635	0.950000000000000	0	0	0.950000000000000
            3636	0.950000000000000	0	0	0.950000000000000
            3637	0.950000000000000	0	0	0.950000000000000
            3638	0.950000000000000	0	0	0.950000000000000
            3639	0.950000000000000	0	0	0.950000000000000
            3640	0.950000000000000	0	0	0.950000000000000
            3641	0.950000000000000	0	0	0.950000000000000
            3642	0.950000000000000	0	0	0.950000000000000
            3643	0.950000000000000	0	0	0.950000000000000
            3644	0.950000000000000	0	0	0.950000000000000
            3645	0.758686982128141	0	0	0.758686982128141
            3646	0.942201652382293	0	0	0.942201652382293
            3647	0.451158294027635	0	0	0.451158294027635
            3648	0.364077142963561	0	0	0.364077142963561
            3649	0.145743466013042	0	0	0.145743466013042
            3650	0.179498279055307	0	0	0.179498279055307
            3651	0.0593369511971704	0	0	0.0593369511971704
            3652	0.0500000000000000	0	0	0.0500000000000000
            3653	0.0500000000000000	0	0	0.0500000000000000
            3654	0.0500000000000000	0	0	0.0500000000000000
            3655	0.0500000000000000	0	0	0.0500000000000000
            3656	0.0500000000000000	0	0	0.0500000000000000
            3657	0.0500000000000000	0	0	0.0500000000000000
            3658	0.0500000000000000	0	0	0.0500000000000000
            3659	0.0500000000000000	0	0	0.0500000000000000
            3660	0.0500000000000000	0	0	0.0500000000000000
            3661	0.0500000000000000	0	0	0.0500000000000000
            3662	0.0500000000000000	0	0	0.0500000000000000
            3663	0.0500000000000000	0	0	0.0500000000000000
            3664	0.0500000000000000	0	0	0.0500000000000000
            3665	0.0500000000000000	0	0	0.0500000000000000
            3666	0.0500000000000000	0	0	0.0500000000000000
            3667	0.0789427048487616	0	0	0.0789427048487616
            3668	0.0870561078190915	0	0	0.0870561078190915
            3669	0.209618199727865	0	0	0.209618199727865
            3670	0.190110156792446	0	0	0.190110156792446
            3671	0.376603734584098	0	0	0.376603734584098
            3672	0.415239028515315	0	0	0.415239028515315
            3673	0.676634437885715	0	0	0.676634437885715
            3674	0.613680327468730	0	0	0.613680327468730
            3675	0.822577580548883	0	0	0.822577580548883
            3676	0.906960750718653	0	0	0.906960750718653
            3677	0.950000000000000	0	0	0.950000000000000
            3678	0.906960637302860	0	0	0.906960637302860
            3679	0.822577749981734	0	0	0.822577749981734
            3680	0.906960757625588	0	0	0.906960757625588
            3681	0.676635132308598	0	0	0.676635132308598
            3682	0.613682020196865	0	0	0.613682020196865
            3683	0.376617209909009	0	0	0.376617209909009
            3684	0.415246825673092	0	0	0.415246825673092
            3685	0.209675334500414	0	0	0.209675334500414
            3686	0.190197330019889	0	0	0.190197330019889
            3687	0.0794041139536000	0	0	0.0794041139536000
            3688	0.0873748941946479	0	0	0.0873748941946479
            3689	0.0500000000000000	0	0	0.0500000000000000
            3690	0.0500000000000000	0	0	0.0500000000000000
            3691	0.0500000000000000	0	0	0.0500000000000000
            3692	0.0500000000000000	0	0	0.0500000000000000
            3693	0.0500000000000000	0	0	0.0500000000000000
            3694	0.0500000000000000	0	0	0.0500000000000000
            3695	0.0513851943122364	0	0	0.0513851943122364
            3696	0.0500000000000000	0	0	0.0500000000000000
            3697	0.0730674801174143	0	0	0.0730674801174143
            3698	0.102022667047169	0	0	0.102022667047169
            3699	0.168184231174880	0	0	0.168184231174880
            3700	0.119544442552393	0	0	0.119544442552393
            3701	0.161349926854857	0	0	0.161349926854857
            3702	0.229417875930213	0	0	0.229417875930213
            3703	0.259456459512407	0	0	0.259456459512407
            3704	0.179728607890625	0	0	0.179728607890625
            3705	0.165456430145824	0	0	0.165456430145824
            3706	0.244035734804927	0	0	0.244035734804927
            3707	0.191675169951791	0	0	0.191675169951791
            3708	0.126139033039976	0	0	0.126139033039976
            3709	0.0798587509974213	0	0	0.0798587509974213
            3710	0.126346261780106	0	0	0.126346261780106
            3711	0.0702843099949103	0	0	0.0702843099949103
            3712	0.0500000000000000	0	0	0.0500000000000000
            3713	0.209611915710572	0	0	0.209611915710572
            3714	0.415237153047482	0	0	0.415237153047482
            3715	0.376606516289529	0	0	0.376606516289529
            3716	0.190114268608497	0	0	0.190114268608497
            3717	0.141868994997792	0	0	0.141868994997792
            3718	0.280987125289903	0	0	0.280987125289903
            3719	0.172574042892395	0	0	0.172574042892395
            3720	0.0872867364091474	0	0	0.0872867364091474
            3721	0.0500000000000000	0	0	0.0500000000000000
            3722	0.0878402876714708	0	0	0.0878402876714708
            3723	0.0500000000000000	0	0	0.0500000000000000
            3724	0.0500000000000000	0	0	0.0500000000000000
            3725	0.0500000000000000	0	0	0.0500000000000000
            3726	0.0500000000000000	0	0	0.0500000000000000
            3727	0.0500000000000000	0	0	0.0500000000000000
            3728	0.0608922740965812	0	0	0.0608922740965812
            3729	0.143541014864516	0	0	0.143541014864516
            3730	0.0886415589635828	0	0	0.0886415589635828
            3731	0.176042412545325	0	0	0.176042412545325
            3732	0.286692018988289	0	0	0.286692018988289
            3733	0.477961185789239	0	0	0.477961185789239
            3734	0.293331789070080	0	0	0.293331789070080
            3735	0.412748056561566	0	0	0.412748056561566
            3736	0.672574907690689	0	0	0.672574907690689
            3737	0.818464219924782	0	0	0.818464219924782
            3738	0.502275586136177	0	0	0.502275586136177
            3739	0.549025865471550	0	0	0.549025865471550
            3740	0.894644156953760	0	0	0.894644156953760
            3741	0.915672551571991	0	0	0.915672551571991
            3742	0.561933925066364	0	0	0.561933925066364
            3743	0.549083430907720	0	0	0.549083430907720
            3744	0.894690488860317	0	0	0.894690488860317
            3745	0.818777770442584	0	0	0.818777770442584
            3746	0.502659258647015	0	0	0.502659258647015
            3747	0.414840378601287	0	0	0.414840378601287
            3748	0.674316521865316	0	0	0.674316521865316
            3749	0.486017505470061	0	0	0.486017505470061
            3750	0.302805166023792	0	0	0.302805166023792
            3751	0.211741714693750	0	0	0.211741714693750
            3752	0.317808255408344	0	0	0.317808255408344
            3753	0.244144058564133	0	0	0.244144058564133
            3754	0.200900747291949	0	0	0.200900747291949
            3755	0.335141315674396	0	0	0.335141315674396
            3756	0.333838416369746	0	0	0.333838416369746
            3757	0.649107939027049	0	0	0.649107939027049
            3758	0.676733579694532	0	0	0.676733579694532
            3759	0.950000000000000	0	0	0.950000000000000
            3760	0.950000000000000	0	0	0.950000000000000
            3761	0.950000000000000	0	0	0.950000000000000
            3762	0.950000000000000	0	0	0.950000000000000
            3763	0.950000000000000	0	0	0.950000000000000
            3764	0.950000000000000	0	0	0.950000000000000
            3765	0.950000000000000	0	0	0.950000000000000
            3766	0.950000000000000	0	0	0.950000000000000
            3767	0.950000000000000	0	0	0.950000000000000
            3768	0.950000000000000	0	0	0.950000000000000
            3769	0.950000000000000	0	0	0.950000000000000
            3770	0.950000000000000	0	0	0.950000000000000
            3771	0.950000000000000	0	0	0.950000000000000
            3772	0.950000000000000	0	0	0.950000000000000
            3773	0.950000000000000	0	0	0.950000000000000
            3774	0.916975651768264	0	0	0.916975651768264
            3775	0.437304583936388	0	0	0.437304583936388
            3776	0.479326855030525	0	0	0.479326855030525
            3777	0.190553397523524	0	0	0.190553397523524
            3778	0.173630349237126	0	0	0.173630349237126
            3779	0.0572544223085148	0	0	0.0572544223085148
            3780	0.0628609303579899	0	0	0.0628609303579899
            3781	0.0500000000000000	0	0	0.0500000000000000
            3782	0.0500000000000000	0	0	0.0500000000000000
            3783	0.0500000000000000	0	0	0.0500000000000000
            3784	0.0500000000000000	0	0	0.0500000000000000
            3785	0.0500000000000000	0	0	0.0500000000000000
            3786	0.0500000000000000	0	0	0.0500000000000000
            3787	0.0500000000000000	0	0	0.0500000000000000
            3788	0.0500000000000000	0	0	0.0500000000000000
            3789	0.0500000000000000	0	0	0.0500000000000000
            3790	0.0500000000000000	0	0	0.0500000000000000
            3791	0.0500000000000000	0	0	0.0500000000000000
            3792	0.0500000000000000	0	0	0.0500000000000000
            3793	0.0500000000000000	0	0	0.0500000000000000
            3794	0.0500000000000000	0	0	0.0500000000000000
            3795	0.0500000000000000	0	0	0.0500000000000000
            3796	0.0588933673591777	0	0	0.0588933673591777
            3797	0.141830254125662	0	0	0.141830254125662
            3798	0.0870383765085163	0	0	0.0870383765085163
            3799	0.172421626568489	0	0	0.172421626568489
            3800	0.280963322659908	0	0	0.280963322659908
            3801	0.457833369969342	0	0	0.457833369969342
            3802	0.280963293356902	0	0	0.280963293356902
            3803	0.376603454170875	0	0	0.376603454170875
            3804	0.613680255057898	0	0	0.613680255057898
            3805	0.676633875316472	0	0	0.676633875316472
            3806	0.415236874903802	0	0	0.415236874903802
            3807	0.376603981793247	0	0	0.376603981793247
            3808	0.613680582557711	0	0	0.613680582557711
            3809	0.457836494819794	0	0	0.457836494819794
            3810	0.280968316170027	0	0	0.280968316170027
            3811	0.172460731038530	0	0	0.172460731038530
            3812	0.280987670634688	0	0	0.280987670634688
            3813	0.141986257870466	0	0	0.141986257870466
            3814	0.0872889463195492	0	0	0.0872889463195492
            3815	0.0500000000000000	0	0	0.0500000000000000
            3816	0.0597158210661209	0	0	0.0597158210661209
            3817	0.0500000000000000	0	0	0.0500000000000000
            3818	0.0500000000000000	0	0	0.0500000000000000
            3819	0.0500000000000000	0	0	0.0500000000000000
            3820	0.0500000000000000	0	0	0.0500000000000000
            3821	0.0500000000000000	0	0	0.0500000000000000
            3822	0.0616135705764241	0	0	0.0616135705764241
            3823	0.148618777307702	0	0	0.148618777307702
            3824	0.0916251187490541	0	0	0.0916251187490541
            3825	0.183504531831360	0	0	0.183504531831360
            3826	0.300032951900433	0	0	0.300032951900433
            3827	0.503079571460393	0	0	0.503079571460393
            3828	0.304864320374656	0	0	0.304864320374656
            3829	0.420321192475169	0	0	0.420321192475169
            3830	0.702469789089330	0	0	0.702469789089330
            3831	0.820394487613562	0	0	0.820394487613562
            3832	0.482482798684911	0	0	0.482482798684911
            3833	0.463119503166537	0	0	0.463119503166537
            3834	0.805540372298373	0	0	0.805540372298373
            3835	0.668750959644071	0	0	0.668750959644071
            3836	0.373653460034568	0	0	0.373653460034568
            3837	0.254830393486690	0	0	0.254830393486690
            3838	0.471874161320426	0	0	0.471874161320426
            3839	0.284064963340538	0	0	0.284064963340538
            3840	0.147682319003897	0	0	0.147682319003897
            3841	0.676634009900954	0	0	0.676634009900954
            3842	0.906960685878052	0	0	0.906960685878052
            3843	0.822578204979047	0	0	0.822578204979047
            3844	0.613681798700990	0	0	0.613681798700990
            3845	0.457845392575362	0	0	0.457845392575362
            3846	0.613685246831824	0	0	0.613685246831824
            3847	0.376635399469659	0	0	0.376635399469659
            3848	0.281040233509718	0	0	0.281040233509718
            3849	0.142234967994794	0	0	0.142234967994794
            3850	0.190277364774496	0	0	0.190277364774496
            3851	0.0796680299547251	0	0	0.0796680299547251
            3852	0.0606454987910325	0	0	0.0606454987910325
            3853	0.0500000000000000	0	0	0.0500000000000000
            3854	0.0500000000000000	0	0	0.0500000000000000
            3855	0.0500000000000000	0	0	0.0500000000000000
            3856	0.0500000000000000	0	0	0.0500000000000000
            3857	0.0500000000000000	0	0	0.0500000000000000
            3858	0.0500000000000000	0	0	0.0500000000000000
            3859	0.0500000000000000	0	0	0.0500000000000000
            3860	0.0890375362950624	0	0	0.0890375362950624
            3861	0.148100574282628	0	0	0.148100574282628
            3862	0.0615327394752636	0	0	0.0615327394752636
            3863	0.0865227492227193	0	0	0.0865227492227193
            3864	0.208358506399162	0	0	0.208358506399162
            3865	0.253548919714387	0	0	0.253548919714387
            3866	0.105283372044468	0	0	0.105283372044468
            3867	0.115082696168129	0	0	0.115082696168129
            3868	0.277148502682143	0	0	0.277148502682143
            3869	0.283668118164928	0	0	0.283668118164928
            3870	0.117793134931619	0	0	0.117793134931619
            3871	0.115136198736256	0	0	0.115136198736256
            3872	0.277209130238093	0	0	0.277209130238093
            3873	0.253948909869841	0	0	0.253948909869841
            3874	0.105634025121078	0	0	0.105634025121078
            3875	0.0883946842301024	0	0	0.0883946842301024
            3876	0.210512233175629	0	0	0.210512233175629
            3877	0.157701398348282	0	0	0.157701398348282
            3878	0.0697862600982932	0	0	0.0697862600982932
            3879	0.0672934815640846	0	0	0.0672934815640846
            3880	0.124539269320864	0	0	0.124539269320864
            3881	0.154768687298680	0	0	0.154768687298680
            3882	0.110986856148906	0	0	0.110986856148906
            3883	0.242308502784867	0	0	0.242308502784867
            3884	0.303381976691365	0	0	0.303381976691365
            3885	0.623600715923103	0	0	0.623600715923103
            3886	0.500270300464646	0	0	0.500270300464646
            3887	0.887213369833854	0	0	0.887213369833854
            3888	0.950000000000000	0	0	0.950000000000000
            3889	0.950000000000000	0	0	0.950000000000000
            3890	0.950000000000000	0	0	0.950000000000000
            3891	0.950000000000000	0	0	0.950000000000000
            3892	0.950000000000000	0	0	0.950000000000000
            3893	0.950000000000000	0	0	0.950000000000000
            3894	0.950000000000000	0	0	0.950000000000000
            3895	0.950000000000000	0	0	0.950000000000000
            3896	0.950000000000000	0	0	0.950000000000000
            3897	0.950000000000000	0	0	0.950000000000000
            3898	0.950000000000000	0	0	0.950000000000000
            3899	0.876811195876919	0	0	0.876811195876919
            3900	0.950000000000000	0	0	0.950000000000000
            3901	0.723102339682048	0	0	0.723102339682048
            3902	0.490879011891200	0	0	0.490879011891200
            3903	0.230974615903518	0	0	0.230974615903518
            3904	0.343019102171868	0	0	0.343019102171868
            3905	0.135780175136108	0	0	0.135780175136108
            3906	0.0909446666816881	0	0	0.0909446666816881
            3907	0.0500000000000000	0	0	0.0500000000000000
            3908	0.0500000000000000	0	0	0.0500000000000000
            3909	0.0500000000000000	0	0	0.0500000000000000
            3910	0.0500000000000000	0	0	0.0500000000000000
            3911	0.0500000000000000	0	0	0.0500000000000000
            3912	0.0500000000000000	0	0	0.0500000000000000
            3913	0.0500000000000000	0	0	0.0500000000000000
            3914	0.0500000000000000	0	0	0.0500000000000000
            3915	0.0500000000000000	0	0	0.0500000000000000
            3916	0.0500000000000000	0	0	0.0500000000000000
            3917	0.0500000000000000	0	0	0.0500000000000000
            3918	0.0500000000000000	0	0	0.0500000000000000
            3919	0.0500000000000000	0	0	0.0500000000000000
            3920	0.0500000000000000	0	0	0.0500000000000000
            3921	0.0500000000000000	0	0	0.0500000000000000
            3922	0.0500000000000000	0	0	0.0500000000000000
            3923	0.0500000000000000	0	0	0.0500000000000000
            3924	0.0500000000000000	0	0	0.0500000000000000
            3925	0.0500000000000000	0	0	0.0500000000000000
            3926	0.0500000000000000	0	0	0.0500000000000000
            3927	0.0500000000000000	0	0	0.0500000000000000
            3928	0.0870383678770271	0	0	0.0870383678770271
            3929	0.141830159396993	0	0	0.141830159396993
            3930	0.0588931057601520	0	0	0.0588931057601520
            3931	0.0789403765860326	0	0	0.0789403765860326
            3932	0.190109277638120	0	0	0.190109277638120
            3933	0.209611448959024	0	0	0.209611448959024
            3934	0.0870384358426311	0	0	0.0870384358426311
            3935	0.0789411561955846	0	0	0.0789411561955846
            3936	0.190109983984409	0	0	0.190109983984409
            3937	0.141836885610143	0	0	0.141836885610143
            3938	0.0589005343385997	0	0	0.0589005343385997
            3939	0.0500000000000000	0	0	0.0500000000000000
            3940	0.0870907643552231	0	0	0.0870907643552231
            3941	0.0500000000000000	0	0	0.0500000000000000
            3942	0.0500000000000000	0	0	0.0500000000000000
            3943	0.0500000000000000	0	0	0.0500000000000000
            3944	0.0500000000000000	0	0	0.0500000000000000
            3945	0.0500000000000000	0	0	0.0500000000000000
            3946	0.0500000000000000	0	0	0.0500000000000000
            3947	0.0500000000000000	0	0	0.0500000000000000
            3948	0.0500000000000000	0	0	0.0500000000000000
            3949	0.0828453575896466	0	0	0.0828453575896466
            3950	0.0928631861612792	0	0	0.0928631861612792
            3951	0.229771153799570	0	0	0.229771153799570
            3952	0.202607792334262	0	0	0.202607792334262
            3953	0.412909407018194	0	0	0.412909407018194
            3954	0.474158760169732	0	0	0.474158760169732
            3955	0.818503480419388	0	0	0.818503480419388
            3956	0.700968489044262	0	0	0.700968489044262
            3957	0.950000000000000	0	0	0.950000000000000
            3958	0.950000000000000	0	0	0.950000000000000
            3959	0.950000000000000	0	0	0.950000000000000
            3960	0.950000000000000	0	0	0.950000000000000
            3961	0.950000000000000	0	0	0.950000000000000
            3962	0.950000000000000	0	0	0.950000000000000
            3963	0.950000000000000	0	0	0.950000000000000
            3964	0.950000000000000	0	0	0.950000000000000
            3965	0.753359292186743	0	0	0.753359292186743
            3966	0.950000000000000	0	0	0.950000000000000
            3967	0.667456778702364	0	0	0.667456778702364
            3968	0.471149492460556	0	0	0.471149492460556
            3969	0.950000000000000	0	0	0.950000000000000
            3970	0.906960624411827	0	0	0.906960624411827
            3971	0.822577624065134	0	0	0.822577624065134
            3972	0.906960837370409	0	0	0.906960837370409
            3973	0.676635552493578	0	0	0.676635552493578
            3974	0.613680730615319	0	0	0.613680730615319
            3975	0.376606516741200	0	0	0.376606516741200
            3976	0.415247741244654	0	0	0.415247741244654
            3977	0.209668801263990	0	0	0.209668801263990
            3978	0.190125404465705	0	0	0.190125404465705
            3979	0.0790102029514600	0	0	0.0790102029514600
            3980	0.0872869098803079	0	0	0.0872869098803079
            3981	0.0500000000000000	0	0	0.0500000000000000
            3982	0.0500000000000000	0	0	0.0500000000000000
            3983	0.0500000000000000	0	0	0.0500000000000000
            3984	0.0500000000000000	0	0	0.0500000000000000
            3985	0.0500000000000000	0	0	0.0500000000000000
            3986	0.0500000000000000	0	0	0.0500000000000000
            3987	0.0500000000000000	0	0	0.0500000000000000
            3988	0.0500000000000000	0	0	0.0500000000000000
            3989	0.0500000000000000	0	0	0.0500000000000000
            3990	0.0500000000000000	0	0	0.0500000000000000
            3991	0.0500000000000000	0	0	0.0500000000000000
            3992	0.0500000000000000	0	0	0.0500000000000000
            3993	0.0500000000000000	0	0	0.0500000000000000
            3994	0.0500000000000000	0	0	0.0500000000000000
            3995	0.0500000000000000	0	0	0.0500000000000000
            3996	0.0500000000000000	0	0	0.0500000000000000
            3997	0.0500000000000000	0	0	0.0500000000000000
            3998	0.0500000000000000	0	0	0.0500000000000000
            3999	0.0500000000000000	0	0	0.0500000000000000
            4000	0.0500000000000000	0	0	0.0500000000000000
            4001	0.0500000000000000	0	0	0.0500000000000000
            4002	0.0500000000000000	0	0	0.0500000000000000
            4003	0.0500000000000000	0	0	0.0500000000000000
            4004	0.0500000000000000	0	0	0.0500000000000000
            4005	0.0500000000000000	0	0	0.0500000000000000
            4006	0.0500000000000000	0	0	0.0500000000000000
            4007	0.0500000000000000	0	0	0.0500000000000000
            4008	0.0500000000000000	0	0	0.0500000000000000
            4009	0.0719650746726178	0	0	0.0719650746726178
            4010	0.0500000000000000	0	0	0.0500000000000000
            4011	0.0974714398647733	0	0	0.0974714398647733
            4012	0.166936577535217	0	0	0.166936577535217
            4013	0.343810271134886	0	0	0.343810271134886
            4014	0.199733823748654	0	0	0.199733823748654
            4015	0.343546305138056	0	0	0.343546305138056
            4016	0.599506040389421	0	0	0.599506040389421
            4017	0.881193930438502	0	0	0.881193930438502
            4018	0.496118409179154	0	0	0.496118409179154
            4019	0.604149076285006	0	0	0.604149076285006
            4020	0.950000000000000	0	0	0.950000000000000
            4021	0.950000000000000	0	0	0.950000000000000
            4022	0.623405648186586	0	0	0.623405648186586
            4023	0.547267186147243	0	0	0.547267186147243
            4024	0.950000000000000	0	0	0.950000000000000
            4025	0.795400244592619	0	0	0.795400244592619
            4026	0.409609926173675	0	0	0.409609926173675
            4027	0.261328986988605	0	0	0.261328986988605
            4028	0.517828096819983	0	0	0.517828096819983
            4029	0.285700501544606	0	0	0.285700501544606
            4030	0.141756151467816	0	0	0.141756151467816
            4031	0.0651007641718449	0	0	0.0651007641718449
            4032	0.132974076541992	0	0	0.132974076541992
            4033	0.0519566004073874	0	0	0.0519566004073874
            4034	0.0500000000000000	0	0	0.0500000000000000
            4035	0.0500000000000000	0	0	0.0500000000000000
            4036	0.0500000000000000	0	0	0.0500000000000000
            4037	0.0500000000000000	0	0	0.0500000000000000
            4038	0.0500000000000000	0	0	0.0500000000000000
            4039	0.0500000000000000	0	0	0.0500000000000000
            4040	0.0500000000000000	0	0	0.0500000000000000
            4041	0.0500000000000000	0	0	0.0500000000000000
            4042	0.0500000000000000	0	0	0.0500000000000000
            4043	0.0500000000000000	0	0	0.0500000000000000
            4044	0.0500000000000000	0	0	0.0500000000000000
            4045	0.0500000000000000	0	0	0.0500000000000000
            4046	0.0500000000000000	0	0	0.0500000000000000
            4047	0.0500000000000000	0	0	0.0500000000000000
            4048	0.0500000000000000	0	0	0.0500000000000000
            4049	0.0500000000000000	0	0	0.0500000000000000
            4050	0.0500000000000000	0	0	0.0500000000000000
            4051	0.0500000000000000	0	0	0.0500000000000000
            4052	0.0500000000000000	0	0	0.0500000000000000
            4053	0.0500000000000000	0	0	0.0500000000000000
            4054	0.0500000000000000	0	0	0.0500000000000000
            4055	0.0500000000000000	0	0	0.0500000000000000
            4056	0.0500000000000000	0	0	0.0500000000000000
            4057	0.0500000000000000	0	0	0.0500000000000000
            4058	0.0500000000000000	0	0	0.0500000000000000
            4059	0.0500000000000000	0	0	0.0500000000000000
            4060	0.0500000000000000	0	0	0.0500000000000000
            4061	0.0500000000000000	0	0	0.0500000000000000
            4062	0.0500000000000000	0	0	0.0500000000000000
            4063	0.0500000000000000	0	0	0.0500000000000000
            4064	0.0500000000000000	0	0	0.0500000000000000
            4065	0.0500000000000000	0	0	0.0500000000000000
            4066	0.0500000000000000	0	0	0.0500000000000000
            4067	0.0500000000000000	0	0	0.0500000000000000
            4068	0.0500000000000000	0	0	0.0500000000000000
            4069	0.0500000000000000	0	0	0.0500000000000000
            4070	0.0500000000000000	0	0	0.0500000000000000
            4071	0.0500000000000000	0	0	0.0500000000000000
            4072	0.0500000000000000	0	0	0.0500000000000000
            4073	0.0500000000000000	0	0	0.0500000000000000
            4074	0.0500000000000000	0	0	0.0500000000000000
            4075	0.0500000000000000	0	0	0.0500000000000000
            4076	0.0500000000000000	0	0	0.0500000000000000
            4077	0.0865740825967350	0	0	0.0865740825967350
            4078	0.0672621167107953	0	0	0.0672621167107953
            4079	0.171562604670348	0	0	0.171562604670348
            4080	0.217084959651261	0	0	0.217084959651261
            4081	0.455548594064585	0	0	0.455548594064585
            4082	0.368002463849415	0	0	0.368002463849415
            4083	0.667445339134614	0	0	0.667445339134614
            4084	0.803794739085117	0	0	0.803794739085117
            4085	0.950000000000000	0	0	0.950000000000000
            4086	0.950000000000000	0	0	0.950000000000000
            4087	0.950000000000000	0	0	0.950000000000000
            4088	0.950000000000000	0	0	0.950000000000000
            4089	0.950000000000000	0	0	0.950000000000000
            4090	0.950000000000000	0	0	0.950000000000000
            4091	0.950000000000000	0	0	0.950000000000000
            4092	0.950000000000000	0	0	0.950000000000000
            4093	0.950000000000000	0	0	0.950000000000000
            4094	0.950000000000000	0	0	0.950000000000000
            4095	0.818464168983747	0	0	0.818464168983747
            4096	0.803791071845297	0	0	0.803791071845297];
        for i = 1:size(centelem,1)
            elem(i,5)=i;
        end
end
kmap=K;
fonte=fonte;
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

function dtetadx=calculotetadx(x,y)

if x>0
    m=(-y)/(x^2+y^2);
elseif x<0 && y>=0
    m=(-y)/(x^2+y^2);
elseif x<0 && y<0
    m=(-y)/(x^2+y^2);
elseif x==0 && y>0
    m=0;
elseif x==0 && y<0
    m=0;
else
    m='indefined';
end
dtetadx=m;
end

function dtetady=calculotetady(x,y)

if x>0
    m=x/(x^2+y^2);
elseif x<0 && y>=0
    m=x/(x^2+y^2);
elseif x<0 && y<0
    m=x/(x^2+y^2);
elseif x==0 && y>0
    m=0;
elseif x==0 && y<0
    m=0;
else
    m='indefined';
end
dtetady=m;
end