function [ lambda,r ] =  Lamdas_Weights_LPEW1( Kt1, Kt2, Kn1, Kn2, theta1,...
    theta2, fi1, fi2, netas, P, O,Qo,T,No,r )


nec=size(O,1);
lambda=zeros(nec,1);

%Cálculo dos psi.%

if size(P,1)==size(O,1) %Se for um nó interno.%%%%%%
    for k=1:nec,
        if (k==1)&&(size(P,1)==size(O,1))
            psin=Kn2(nec,2)*cot(fi2(nec))+Kn2(k,1)*cot(fi1(k))+Kt2(nec,2) ...
                -Kt2(k,1);
            psid=Kn1(nec,2)*cot(theta2(nec))+Kn1(k,1)*cot(theta1(k)) ...
                -Kt1(nec,2)+Kt1(k,1);
        else
            psin=Kn2(k-1,2)*cot(fi2(k-1))+Kn2(k,1)*cot(fi1(k))+Kt2(k-1,2) ...
                -Kt2(k,1);
            psid=Kn1(k-1,2)*cot(theta2(k-1))+Kn1(k,1)*cot(theta1(k))- ...
                Kt1(k-1,2)+Kt1(k,1);
        end
        psi(k)=psin/psid;
    end
    
else %Se for um nó do contorno.%%%%%%%%%
    for k=1:nec+1,
        if (k==1)&&(size(P,1)~=size(O,1))
            psin=Kn2(k,1)*cot(fi1(k))-Kt2(k,1);
            psid=Kn1(k,1)*cot(theta1(k))+Kt1(k,1);
            % calculo do r
            r(No,1)=(1+ (psin/psid))*norm(Qo-T(1,:));
        elseif (k==nec+1)&&(size(P,1)~=size(O,1))
            psin=Kn2(k-1,2)*cot(fi2(k-1))+Kt2(k-1,2);
            psid=Kn1(k-1,2)*cot(theta2(k-1))-Kt1(k-1,2);
            % calculo do r
            
            r(No,2)=(1+(psin/psid))*norm(Qo-T(nec+1,:));
        else
            psin=Kn2(k-1,2)*cot(fi2(k-1))+Kn2(k,1)*cot(fi1(k))+Kt2(k-1,2) ...
                -Kt2(k,1);
            psid=Kn1(k-1,2)*cot(theta2(k-1))+Kn1(k,1)*cot(theta1(k))- ...
                Kt1(k-1,2)+Kt1(k,1);
        end
        psi(k)=psin/psid;
    end
end

%Cálculo dos lambdas.%

for k=1:nec,
    if (k==nec)&&(size(P,1)==size(O,1))
        lambda(k)=Kn1(k,1)*netas(k,1)*psi(k)+Kn1(k,2)*netas(k,2)* ...
            psi(1)-(Kn2(k,1)*cot(theta1(k)+fi1(k))+Kn2(k,2)* ...
            cot(theta2(k)+fi2(k)))+Kt2(k,1)-Kt2(k,2);
    else
        lambda(k)=Kn1(k,1)*netas(k,1)*psi(k)+Kn1(k,2)*netas(k,2)* ...
            psi(k+1)-(Kn2(k,1)*cot(theta1(k)+fi1(k))+Kn2(k,2)* ...
            cot(theta2(k)+fi2(k)))+Kt2(k,1)-Kt2(k,2);
    end
end
end

