function [ w,ww] = Weights_Interp_LPEW2mod(kmap)

global esurn2 coord

%Retorna os valores dos pesos w'e w%

%Prealocação dos vetores%
%ww=zeros(esurn2(ni+1)-esurn2(ni),2);
%w=zeros(esurn2(ni+1)-esurn2(ni),2);

%Equação Base%
%w'(k)=(((netas(k,1)+netas(k-1,1))*E(k,1))/(Ë(k,1)+Ë(k-1,1)))+(((neta(k,2)+neta(k+1,2))*E(k,2))/(Ë(k,2)+Ë(k+1,2)))%

%Loop que percorre os elementos em torno do nó "ni".%
for No=1:size(coord,1)
    
    [ O, P, T, Qo, areas, distancias, talfa] = OPTADTa_Interp_LPEW2mod(No);
    [ netas ] =netas_Interp_LPEW2mod(O, P, T, Qo, No, kmap);
   
    [ E,EE ] =EE_Interp_LPEW2mod(O, P, T, Qo, No, kmap);
    
    
    for k=1:size(O,1)
        
        ww(k)=(((netas(k,1)+netas(k-1,1))*E(k,1))/(EE(k,1)+EE(k-1,1)))+(((neta(k,2)+neta(k+1,2))*E(k,2))/(EE(k,2)+EE(k+1,2)))
        
    end
    
    for k=0:size(O,1)-1
        
        w(k)=ww(k+1)/sum(ww);
        
    end
end
end




