% Esta funcao calcula os pesos fornecidos na equacao 15do artigo  S.Miao
% and J. Wu, 2021.
function [ w,s] = Pre_ELPEW_2(kmap)
global coord  esurn2
% devolve os pesos "w" cujos elementos são organizados em um vetor
%Retorna todos os parâmetros necessários às expressões dos fluxos.
apw=ones(size(coord,1),1);
r=zeros(size(coord,1),2);
s=0;
epsilon=0.5;

for No=1:size(coord,1)
    % calcula
    
   O=zeros(esurn2(No+1)-esurn2(No),3); 
   
   [netas,E,EE]=parameters_weight(No,kmap,epsilon);
   [omegas]=aux_omegacalculate(E,EE,netas,No);
    for k=0:size(O,1)-1
        w(apw(No)+k,1)=omegas(k+1)/sum(omegas); %Os pesos fazem sentido
    end
    apw(No+1)=apw(No)+size(O,1);
end
end

