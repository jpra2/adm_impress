function[N,F,V,weightLS,esuel,esuel_coord,A,bound,S_old,S_cont]=presaturation(wells)

global satlimit elem bedge bcflag
% adequa��o das faces por elemento e por n� em sentido anti-horario
% F: vetor de faces em cada elemento
% V: faces ordenados "convenientemente" ao rededor de um n�
% N: faces ordenados ao rededor de um no em sentido anti-horario

[F,V,N]=elementface;

[esuel, esuel_coord,A,bound] = esurnelem;

% calculo dos pesos para o m�todo minimos quadrados
[weightLS] = weights(A,esuel,esuel_coord);

% Condi�ao inicial satura�ao
S_old = zeros(size(elem,1),1);

S_old(:)=satlimit(2);
% Condi�ao de contorno da satura�ao

S_cont=1-satlimit(1);

% adequa��o dos po�os
if max(max(wells))~=0
    for i=1:size(wells,1)
        if wells(i,2)==1
            S_old(wells(i,1))=S_cont;
        end
    end
else
    for ifacont=1:size(bedge,1)
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        rr=bcflag(r,2);
        
        if rr~=0 && bcflag(r,1)>200
            S_cont=1-satlimit(1);
        end
        
    end
end
end