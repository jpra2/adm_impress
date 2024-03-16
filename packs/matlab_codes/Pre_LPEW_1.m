function [ w,s] = Pre_LPEW_1(kmap,N)
%Retorna todos os parâmetros necessários às expressões dos fluxos.

global coord nsurn1 nsurn2 bcflag bedge inedge
apw=ones(size(coord,1),1);
r=zeros(size(coord,1),2);
for y=1:size(coord,1),
    No=y;
    % calculos dos vetores O, P, T, Q
    [ O, P, T, Qo ] = OPT_Interp_LPEW (No);
    % calculo dos angulos
    [ fi2, fi1, theta2, theta1 ] = angulos_Interp_LPEW1( O, P, T, Qo, No );
    % calculo dos netas
    [ neta ] = netas_Interp_LPEW( O, P, T, Qo,No );
    % calculo dos Ks
    [ Kt1, Kt2, Kn1, Kn2 ] = Ks_Interp_LPEW1( O, T, Qo, kmap, No);
    % calculo dos lamdas
    [ lambda,r] = Lamdas_Weights_LPEW1( Kt1, Kt2, Kn1, Kn2, theta1, ...
        theta2, fi1, fi2, neta, P, O,Qo,T,No,r);
    % calculo dos pesos
    for k=0:size(O,1)-1,
        w(apw(No)+k)=lambda(k+1)/sum(lambda); 
    end
    apw(No+1)=apw(No)+size(O,1);
    % interpolaçao das pressões nos contornos de Neumann
    vetor=nsurn1(nsurn2(No)+1:nsurn2(No+1));
    comp1=N(No,1);
    comp2=N(No,length(vetor));
    if comp1>size(inedge,1) && comp2>size(inedge,1)
        a=bcflag(:,1)==bedge(comp1-size(inedge,1),5);
        s1=find(a==1);
        b=bcflag(:,1)==bedge(comp2-size(inedge,1),5);
        s2=find(b==1);
        
        s(No,1)=-(1/sum(lambda))*(r(No,1)*bcflag(s1,2)+r(No,2)*bcflag(s2,2));
    end
end


end

