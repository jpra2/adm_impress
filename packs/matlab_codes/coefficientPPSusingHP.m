function [coefficient]=coefficientPPSusingHP(kmap,facelement,harmonicpoint)
global inedge bedge coord elem centelem
Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0;0 0 0];


for ifacont=1:size(bedge,1)
    
    lef=bedge(ifacont,3);
    % usando quando houver malha malhas mistas
    klef=elementype(lef);
    
    IJ=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    normIJ=norm(IJ);
    %Essa é UMA maneira de construir os tensores.
    
    Klef(1,1)=kmap(elem(bedge(ifacont,3),5),2);
    Klef(1,2)=kmap(elem(bedge(ifacont,3),5),3);
    Klef(2,1)=kmap(elem(bedge(ifacont,3),5),4);
    Klef(2,2)=kmap(elem(bedge(ifacont,3),5),5);
    
    % ej(iface) do elemento a esquerda
    ve2= Klef'*(R*(IJ/normIJ)');
    % os vértices que compõem o elemento lef, ordenados em sentido
    % antihorario
    auxvetor=elem(lef,1:klef);
    
    j=1;
    for i=auxvetor
        if j==length(auxvetor)
            no1=i;
            no2=auxvetor(1);
            % alfa i correspondente ao no1
            % alfa j correspondente ao no2
            [alfai,alfaj]=calconormal(no1,no2,Klef,lef);
            
            aproxerro=norm((coord(no1,:)-centelem(lef,:))*alfai+ (coord(no2,:)-centelem(lef,:))*alfaj-ve2');
            if aproxerro<1e-10 && alfaj>0 && alfai>0
                
                ksii=alfai;
                ksij=alfaj;
                % indexando as faces respetivamente
                aux11=no1;
                aux12=no2;
                % verificando a identidade
                aux13=aproxerro;
            end
        else
            no1=i;
            no2=auxvetor(j+1);
            % alfa i correspondente ao no1
            % alfa j correspondente ao no2
            [alfai,alfaj]=calconormal(no1,no2,Klef,lef);
            
            
            aproxerro=norm((coord(no1,:)-centelem(lef,:))*alfai+ (coord(no2,:)-centelem(lef,:))*alfaj-ve2');
            
            if aproxerro<1e-10 && alfaj>0 && alfai>0
                
                ksii=alfai;
                ksij=alfaj;
                % indexando as faces respetivamente
                aux11=no1;
                aux12=no2;
                % verificando a identidade
                aux13=aproxerro;
                
            end
        end
        j=j+1;
    end
    
    coefficient(1,1,ifacont)=ksii;
    coefficient(1,2,ifacont)=ksij;
    % indexando as faces respetivamente
    coefficient(1,3,ifacont)=aux11;
    coefficient(1,4,ifacont)=aux12;
    % verificando a identidade
    coefficient(1,5,ifacont)=aux13;
    
    clear ksii ksij  aux11 aux12 aux13
end

%% Faces interiores
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    klef=elementype(lef);
    krel=elementype(rel);
    IJ=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    normIJ=norm(coord(inedge(iface,2),:)-coord(inedge(iface,1),:));
    %Essa é UMA maneira de construir os tensores.
    
    Klef(1,1)=kmap(elem(inedge(iface,3),5),2);
    Klef(1,2)=kmap(elem(inedge(iface,3),5),3);
    Klef(2,1)=kmap(elem(inedge(iface,3),5),4);
    Klef(2,2)=kmap(elem(inedge(iface,3),5),5);
    
    %tensor a elemento a direita
    
    Krel(1,1)=kmap(elem(inedge(iface,4),5),2);
    Krel(1,2)=kmap(elem(inedge(iface,4),5),3);
    Krel(2,1)=kmap(elem(inedge(iface,4),5),4);
    Krel(2,2)=kmap(elem(inedge(iface,4),5),5);
    
    %% ej(iface) do elemento a esquerda
    % terceira maneira de calcular o kn
    ve2=Klef*(R*(IJ/normIJ)');
    % os vértices que compõem o elemento lef, ordenados em sentido
    % antihorario
    auxvetor=elem(lef,1:klef);
    
    j=1;
    for i=auxvetor
        if j==length(auxvetor)
            
            no1=i;
            no2=auxvetor(1);
             % alfa i correspondente ao no1
            % alfa j correspondente ao no2
            [alfai,alfaj]=calconormal(no1,no2,Klef,lef);
            
            aproxerro=norm((coord(no1,:)-centelem(lef,:))*alfai+ (coord(no2,:)-centelem(lef,:))*alfaj-ve2');
            if aproxerro<1e-10 && alfaj>0 && alfai>0
                % atribuindo valores a os coeficientes
                ksii=alfai;
                ksij=alfaj;
                % indexando as faces respetivamente
                aux11=no1;
                aux12=no2;
                % verificando a identidade
                aux13=aproxerro;
                
                
            end
            
        else
            no1=i;
            no2=auxvetor(j+1);
            % alfa i correspondente ao no1
            % alfa j correspondente ao no2
            [alfai,alfaj]=calconormal(no1,no2,Klef,lef);
            
            aproxerro=norm((coord(no1,:)-centelem(lef,:))*alfai+ (coord(no2,:)-centelem(lef,:))*alfaj-ve2');
            if aproxerro<1e-10 && alfaj>0 && alfai>0
                
                % atribuindo valores a os coeficientes
                ksii=alfai;
                ksij=alfaj;
                % indexando as faces respetivamente
                aux11=no1;
                aux12=no2;
                % verificando a identidade
                aux13=aproxerro;
            end
        end
        j=j+1;
    end
    % atribuindo valores a os coeficientes
    coefficient(1,1,iface+size(bedge,1))=ksii;
    coefficient(1,2,iface+size(bedge,1))=ksij;
    % indexando as faces respetivamente
    coefficient(1,3,iface+size(bedge,1))=aux11;
    coefficient(1,4,iface+size(bedge,1))=aux12;
    % verificando a identidade
    coefficient(1,5,iface+size(bedge,1))=aux13;
    
    %clear auxi alfaj  no1 no2 aproxerro
    clear ksii ksij  aux11 aux12 aux13
    %% Elemento a direita
    vetor12=Krel*(R*(-IJ/normIJ)');
    % os vértices que compõem o elemento rel, ordenados em sentido
    % antihorario
    auxvetor=elem(rel,1:krel);
    
    j=1;
    for ii=auxvetor
        if j==length(auxvetor)
            no1=ii;
            no2=auxvetor(1);
            
            [alfai,alfaj]=calconormal(no1,no2,Krel,rel);
            
            
            aproxerro=norm((coord(no1,:)-centelem(rel,:))*alfai+ (coord(no2,:)-centelem(rel,:))*alfaj-vetor12');
            if aproxerro<1e-10 && alfaj>0 && alfai>0
                
                ksii=alfai;
                ksij=alfaj;
                aux11=no1;
                aux12=no2;
                aux13=aproxerro;
                
            end
        else
            no1=ii;
            no2=auxvetor(j+1);
            [alfai,alfaj]=calconormal(no1,no2,Krel,rel);
            
            aproxerro=norm((coord(no1,:)-centelem(rel,:))*alfai+ (coord(no2,:)-centelem(rel,:))*alfaj-vetor12');
            if aproxerro<1e-10 && alfaj>0 && alfai>0
                ksii=alfai;
                ksij=alfaj;
                aux11=no1;
                aux12=no2;
                aux13=aproxerro;
                
                
            end
        end
        j=j+1;
        
    end
    % atribuindo valores a os coeficientes
    coefficient(2,1,iface+size(bedge,1))=ksii;
    coefficient(2,2,iface+size(bedge,1))=ksij;
    % indexando as faces respetivamente
    coefficient(2,3,iface+size(bedge,1))=aux11;
    coefficient(2,4,iface+size(bedge,1))=aux12;
    % verificando a identidade
    coefficient(2,5,iface+size(bedge,1))=aux13;
    
    %clear alfai alfaj  no1 no2 aproxerro
    clear ksii ksij  aux11 aux12 aux13
    
end
