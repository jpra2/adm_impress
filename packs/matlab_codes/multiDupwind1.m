function [S_old,RHS,Flux,Q]=multiDupwind1(elem,inedge,bedge,esurn1, ...
    esurn2, bflux, influx,S_old,d_t,pormap, volume,f_elem,w,nsurn2,...
    faces_esurn,f_bedge,coord,vconst1,q,wells,noelement)
RHS=zeros(size(elem,1),1);
Q=zeros(size(elem,1),1);

for iw=1:size(wells,1)
    
    if wells(iw,2)==1 | wells(iw,2)==2
        
        RHS(wells(iw,1)) = RHS(wells(iw,1)) + f_elem(wells(iw,1))*q(wells(iw,1));
    end
end


for i=1:size(inedge,1)
    
    v1= inedge(i,1);
    v2= inedge(i,2);
    lef=inedge(i,3);
    rel=inedge(i,4);
    
    %% Calculo dos fluxo nos vertices
    for ielem =[lef rel]
        
        for ino=[ v1 v2]
            v= bedge(:,1)==ino;
            r=max(v);
            ss=noelement(ino,:)==ielem;
            m=max(ss);
            if r==0 && m~=0
                
                [fluxno,jelem]= fluxnode(esurn1, esurn2,...
                    inedge,ino,influx,w,nsurn2,faces_esurn,ielem);
                
                [fluxno_return,ielemaux]= fluxnode(esurn1, esurn2,...
                    inedge,ino,influx,w,nsurn2,faces_esurn,jelem);
                
                Fluxo= fluxno - fluxno_return;
                
                vmas=max(Fluxo,0);
                vmenos= max(-Fluxo,0);
                
                RHS(ielem)=RHS(ielem) - (f_elem(ielem)*vmas- f_elem(jelem)*vmenos);
                
                RHS(jelem)=RHS(jelem) + (f_elem(ielem)*vmas- f_elem(jelem)*vmenos);
                
                Q(ielem)=Q(ielem) - Fluxo;
                
                Q(jelem)=Q(jelem) + Fluxo;
                
                yy=noelement(ino,:)==jelem;
                noelement(ino,find(yy==1))=0;
                noelement(ino,find(ss==1))=0;
                Zeros=noelement;
            end
        end
    end
    
    

    %% Calculo dos fluxos nas faces interiores
    
    [flux_face]=fluxfaceint(inedge,v1,v2,influx,bedge,w,...
        nsurn2,faces_esurn,i,bflux,lef,rel);
    
    [flux_face_return]=fluxfaceint(inedge,v2,v1,influx,bedge,w,...
        nsurn2,faces_esurn,i,bflux,rel,lef);
    
    Flux(i,1) = flux_face - flux_face_return;
    
    vmas1=max(Flux(i,1),0);
    
    vmenos1=max(-Flux(i,1),0);

    RHS(lef) = RHS(lef)  - (f_elem(lef)*vmas1 - f_elem(rel)*vmenos1);
    
    RHS(rel) = RHS(rel)  + (f_elem(lef)*vmas1 - f_elem(rel)*vmenos1);
    
    Q(lef)= Q(lef)- Flux(i,1);
    Q(rel)= Q(rel)+ Flux(i,1);
end

for i = 1:size(S_old,1)
    porosity=pormap(1);
    if S_old(i,1)~=1
        S_old(i,1) = S_old(i,1) + (d_t*RHS(i))/(porosity*volume(i));
    end
end
end