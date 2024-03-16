function [ lambda,r ] = Lamdas_Weights_LPEW2( Kt1, Kt2, Kn1, Kn2, theta1,...
                                         theta2, ve1, ve2, netas, P, O,Qo,No,T,r )
%Determina os lambdas.
nec=size(O,1);
lambda=zeros(nec,1);

if size(P,1)==size(O,1) %Se for um nó interno.
    for k=1:nec,
        if (k==1)&&(size(P,1)==size(O,1))
            zetan=Kn2(nec)*cot(ve1(nec))+Kn2(k)*cot(ve2(k))+Kt2(nec)-Kt2(k);
            zetad=Kn1(nec,2)*cot(theta2(nec))+Kn1(k,1)*cot(theta1(k)) ...
                -Kt1(nec,2)+Kt1(k,1);
        else
            zetan=Kn2(k-1)*cot(ve1(k-1))+Kn2(k)*cot(ve2(k))+Kt2(k-1)-Kt2(k);
            zetad=Kn1(k-1,2)*cot(theta2(k-1))+Kn1(k,1)*cot(theta1(k)) ...
                -Kt1(k-1,2)+Kt1(k,1);
        end
        zeta(k)=zetan/zetad;    
    end
else %Se for um nó do contorno.
    for k=1:nec+1,
        if (k==1)&&(size(P,1)~=size(O,1))
            zetan=Kn2(k)*cot(ve2(k))-Kt2(k);
            zetad=Kn1(k,1)*cot(theta1(k))+Kt1(k,1);
            r(1,1)=1+ (zetan/zetad);
            % comentei porque ja coloquei a norma em Pre_LPEW2 linha 55
            %r(No,1)=(1+ (zetan/zetad))*norm(Qo-T(1,:));
        elseif (k==nec+1)&&(size(P,1)~=size(O,1))
            zetan=Kn2(k-1)*cot(ve1(k-1))+Kt2(k-1);
            zetad=Kn1(k-1,2)*cot(theta2(k-1))-Kt1(k-1,2);
            r(1,2)=1+(zetan/zetad);
            
            % comentei porque ja coloquei a norma em Pre_LPEW2 linha 55
            %r(No,2)=(1+(zetan/zetad))*norm(Qo-T(nec+1,:));
        else
            zetan=Kn2(k-1)*cot(ve1(k-1))+Kn2(k)*cot(ve2(k))+Kt2(k-1)-Kt2(k);
            zetad=Kn1(k-1,2)*cot(theta2(k-1))+Kn1(k,1)*cot(theta1(k)) ...
                -Kt1(k-1,2)+Kt1(k,1);  
        end
        zeta(k)=zetan/zetad; 
    end
end

for k=1:nec,
    if (k==nec)&&(size(P,1)==size(O,1))
        lambda(k)=Kn1(k,1)*netas(k,1)*zeta(k)+Kn1(k,2)*netas(k,2)*zeta(1);
    else
        lambda(k)=Kn1(k,1)*netas(k,1)*zeta(k)+Kn1(k,2)*netas(k,2)*zeta(k+1);
    end
end

end

