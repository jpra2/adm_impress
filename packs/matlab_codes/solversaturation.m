function [S_old]=solversaturation(S_old,influx,d_t,esuel,wells,totalflow,...
    f_elem,S_cont,W,bound,smetodo,timeorder,upsilon,kappa,nw,no,auxflag)

switch smetodo
    
    case 'FOU'
        [S_old]= firstorderstandard(S_old,influx,totalflow,f_elem,d_t,wells,S_cont,auxflag,nw,no);
        name = smetodo;
        X = sprintf('Calculo do campo de saturação pelo método: %s ',name);
        disp(X)
        
    case {'HOMFV','HOFV-E'}
        
        [Sat_max_min]=Saturation_max_min(S_old);
        switch timeorder
            case 1 % RK!
                [S_old] = highorderstandard(S_old,influx,d_t,esuel,...
                    wells,totalflow,f_elem,S_cont,Sat_max_min,W,bound,upsilon,kappa,nw,no,auxflag,smetodo);
                name = smetodo;
                X = sprintf('Calculo do campo de saturação pelo método: %s ',name);
                disp(X)
            case 2 % RK2
                % paso1
                [S_old1] = highorderstandard(S_old,influx,d_t,esuel,...
                    wells,totalflow,f_elem,S_cont,Sat_max_min,W,bound,upsilon,kappa,nw,no,auxflag,smetodo);
                %paso 2
                [S_old2] = highorderstandard(S_old1,influx,d_t,esuel,...
                    wells,totalflow,f_elem,S_cont,Sat_max_min,W,bound,upsilon,kappa,nw,no,auxflag,smetodo);
                
                S_old=1/2.*(S_old+S_old2);
                disp('Calculo da saturação pelo método MUSCL+RK2')
        end
        
    case 'GOE'
        
        [S_old]=multiDupwind(elem,inedge,bedge,esurn1,esurn2, bflux, influx,...
            S_old,d_t,pormap, volume,f_elem,peso,nsurn2,N,f_bedge,coord,v,...
            wells,totalflux,S_cont,noelemento,elemphant,elembedge,face_bedge);
         
end

end