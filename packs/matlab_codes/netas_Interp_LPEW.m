function [ netas ] =netas_Interp_LPEW(O, P, T, Qo, ni )

global esurn2
% Retorna os netas.
% esta variavel podemos achar no pag. 07 do artigo  chinês

%Prealocação do vetore.%
netas=zeros(esurn2(ni+1)-esurn2(ni),2);

%Loop que percorre os elementos em torno do nó "ni".%

for k=1:size(netas,1),
    
    
    %Preenchimento da segunda coluna do vetor "netas".%
    
    if (k==size(netas,1))&&(size(P,1)==size(O,1))
        v1=O(k,:)-Qo;
        v2=P(1,:)-Qo;
        ce=cross(v1,v2); % produto vetorial
        H2=norm(ce)/norm(v2); % altura
        netas(k,2)=norm(T(1,:)-Qo)/H2;
    else
        v1=O(k,:)-Qo;
        v2=P(k+1,:)-Qo;
        ce=cross(v1,v2);
        H2=norm(ce)/norm(v2); % altura
        netas(k,2)=norm(T(k+1,:)-Qo)/H2;
    end
    
    %%%%%%%Fim do Preenchimento da segunda coluna do vetor "netas".%%%%%%%
    
    
    %Preenchimento da primeira coluna do vetor "netas".%
    
    v1=O(k,:)-Qo;
    v2=P(k,:)-Qo;
    ce=cross(v1,v2); % produto vetorial
    H1=norm(ce)/norm(v2); % altura
    netas(k,1)=norm(T(k,:)-Qo)/H1;
    %%%%%Fim do Preenchimento da primeira coluna do vetor "netas".%%%%%%%

end

