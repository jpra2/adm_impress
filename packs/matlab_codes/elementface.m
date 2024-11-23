function [F,V,N]=elementface
global bedge inedge elem coord esurn1 esurn2 nsurn1 nsurn2

% cuidado quando esta em sentido horario e antihorario.

auxnflag=50000*ones(size(coord,1),1);
for ifacont=1:size(bedge,1)
   
    auxnflag(bedge(ifacont,1),1)=201;
    
end

F=zeros(size(elem,1),4);
%--------------------------------------------------------------------------
% calculate the faces in the neiboring of element  
for ii=1:size(elem,1)
    i=1;
    % decide whether the element is quadrangular or triangular
    list1=4*logical(elem(ii,4)~=0)+3*logical(elem(ii,4)==0);
    % alocate the vertex of element ii
    list=elem(ii,1:list1);
    
    for jj=1:length(list)
        if jj<length(list)
            %ibedge=find(((bedge(:,1)==list(jj+1) & bedge(:,2)==list(jj))|(bedge(:,1)==list(jj) & bedge(:,2)==list(jj+1))));
            b =logical(((bedge(:,1)==list(jj+1) & bedge(:,2)==list(jj))|(bedge(:,1)==list(jj) & bedge(:,2)==list(jj+1))));
            key=1:length(b);
            ibedge=key(b==1);
            if ibedge~=0
                edgeflag=ibedge; % flag da face 
                
            else
                %iedge=find((inedge(:,1)==list(jj+1) & inedge(:,2)==list(jj))|(inedge(:,1)==list(jj) & inedge(:,2)==list(jj+1)));
                b=logical(((inedge(:,1)==list(jj+1) & inedge(:,2)==list(jj))|(inedge(:,1)==list(jj) & inedge(:,2)==list(jj+1))));
                key=1:length(b);
                iedge=key(b==1);
                edgeflag=iedge+size(bedge,1);
                
            end
            
        else
            %ibedge=find(((bedge(:,1)==list(jj) & bedge(:,2)==list(1))|(bedge(:,1)==list(1) & bedge(:,2)==list(jj))));
            b=logical(((bedge(:,1)==list(jj) & bedge(:,2)==list(1))|(bedge(:,1)==list(1) & bedge(:,2)==list(jj))));
            key=1:length(b);
            ibedge=key(b==1);
            if ibedge~=0
                edgeflag=ibedge;
                
            else
                %iedge=find((inedge(:,1)==list(jj) & inedge(:,2)==list(1))|(inedge(:,1)==list(1) & inedge(:,2)==list(jj)));
                
                b=logical((inedge(:,1)==list(jj) & inedge(:,2)==list(1))|(inedge(:,1)==list(1) & inedge(:,2)==list(jj)));
                key=1:length(b);
                iedge=key(b==1);
                edgeflag=iedge+size(bedge,1);
            end
            
        end
        F(ii,i)=edgeflag;
        i=i+1;
    end
    
    
end
for No=1:size(coord,1)
    N_element_No=esurn2(No+1)-esurn2(No);
    
    n_pontos=nsurn2(No+1)-nsurn2(No);
    
    
    
    for k=1:N_element_No
        
        n1=elem(esurn1(esurn2(No)+k),1);
        n2=elem(esurn1(esurn2(No)+k),2);
        n3=elem(esurn1(esurn2(No)+k),3);
        n4=elem(esurn1(esurn2(No)+k),4);
        a=zeros(2,1);
        ii=1;
        for jj=[n1,n2,n3,n4]
            if jj~=No && jj==0
                a(ii,1)=jj;
                ii=ii+1;
            elseif jj~=No && jj~=0
                for g=1:n_pontos
                    h=nsurn1(nsurn2(No)+g);
                    if jj==h
                        a(ii,1)=jj;
                        ii=ii+1;
                    end
                end
            end
        end
        list=nsurn1(nsurn2(No)+1:nsurn2(No+1));
        list2=esurn1(esurn2(No)+1:esurn2(No+1));
        
        for g=1:size(list,1)
           
            h=list(g);
            if length(list)==length(list2)
                if length(list)==k
                    if a(1,1)==list(g)
                        r=size(list,1)+1;
                    elseif a(2,1)==h
                        s=length(list);
                    end
                else
                    if a(1,1)==h
                        r=g;
                    elseif a(2,1)==h
                        s=g;
                    end
                end
            else
                if a(1,1)==h
                    r=g;
                elseif a(2,1)==h
                    s=g;
                end
            end
        end
        % quere dizer que nó pertece ao no interior da malha
        %r=find(bedge(:,1)~=No); 
        if auxnflag(No,1)==50000
            
            if r==n_pontos & s==n_pontos
                
                
                m=a(2,1);
                n=a(1,1);
                
            elseif r<s
                m=a(1,1);
                n=a(2,1);
            else
                m=a(2,1);
                n=a(1,1);
            end
            
            
        else
            if r<s
                m=a(1,1);
                n=a(2,1);
            else
                m=a(2,1);
                n=a(1,1);
            end
        end
        
        
        [Tt]=faces_no(bedge, inedge,No,m,n); % modificar
        V(:,k,No)= Tt';
 
    end
    % N armazena as faces na vizinhanca de um vertice no sentido anti-horario
    m1=1;
    vetor1=nsurn1(nsurn2(No)+1:nsurn2(No+1));
    for j= [vetor1']
        %ibedge=find(((bedge(:,1)==No & bedge(:,2)==j)|(bedge(:,1)==j & bedge(:,2)==No)));
        b=logical(((bedge(:,1)==No & bedge(:,2)==j)|(bedge(:,1)==j & bedge(:,2)==No)));
        key=1:length(b);
        ibedge=key(b==1);
        if ibedge~=0
            
            N(No,m1)=ibedge;
            m1=m1+1;
        else
            %iedge=find(((inedge(:,1)==j & inedge(:,2)==No)|(inedge(:,1)==No & inedge(:,2)==j)));
            b=logical(((inedge(:,1)==j & inedge(:,2)==No)|(inedge(:,1)==No & inedge(:,2)==j)));
            key=1:length(b);
            iedge=key(b==1);
            N(No,m1)=iedge+size(bedge,1);
            m1=m1+1;
        end
    end

end

end