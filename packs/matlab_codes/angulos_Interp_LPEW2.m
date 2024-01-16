function [ ve2, ve1, theta2, theta1 ] =  angulos_Interp_LPEW2( O, P, T, Qo, ni )
global esurn2
%Retorna os �ngulos fi e theta. Ex: fi1=[f1(cell1) fi1(cell2)...
%fi1(celln)], onde n=esurn2(ni+1)-esurn2(ni);

ve2=zeros(1,esurn2(ni+1)-esurn2(ni)); 
ve1=zeros(1,esurn2(ni+1)-esurn2(ni)); 
theta2=zeros(1,esurn2(ni+1)-esurn2(ni));
theta1=zeros(1,esurn2(ni+1)-esurn2(ni)); 

for k=1:size(ve2,2),
%Determina��o dos vetores necess�rios � obten��o dos cossenos:
v0=O(k,:)-Qo;
if (k==size(ve2,2))&&(size(P,1)==size(O,1))
    vetorth2=T(1,:)-Qo;
    vetor1=T(1,:)-T(k,:);
else
    vetor1=T(k+1,:)-T(k,:);
    vetorth2=T(k+1,:)-Qo;
end
vetorth1=T(k,:)-Qo;
%Determina��o dos �ngulos:
ve1(k)=acos(dot(-vetorth1,vetor1)/(norm(vetor1)*norm(vetorth1))); % revisar esses signos
ve2(k)=acos(dot(-vetorth2,-vetor1)/(norm(vetor1)*norm(vetorth2)));
theta2(k)=acos(dot(v0,vetorth2)/(norm(v0)*norm(vetorth2)));
theta1(k)=acos(dot(v0,vetorth1)/(norm(v0)*norm(vetorth1)));
end

end

