function [partialflux]=calfluxopartial2(iface1,iface2,...
              auxparameter1, auxparameter2,gamma,ielem,pinterp,ifactual,p,norma)
y=1000;
if iface1==ifactual
    
    partialflux=norma*(gamma*auxparameter1*(p(ielem)-pinterp(iface1))+ auxparameter2*(p(ielem)-pinterp(iface2)));
    auxparameter=auxparameter1;
elseif iface2==ifactual
    
    partialflux=norma*(gamma*auxparameter2*(p(ielem)-pinterp(iface2))+ auxparameter1*(p(ielem)-pinterp(iface1)));
    auxparameter=auxparameter2;
else
   
    partialflux=norma*(auxparameter1*(p(ielem)-pinterp(iface1))+ auxparameter2*(p(ielem)-pinterp(iface2)));
    auxparameter=0;
    y=0;
end
if abs(partialflux)<1e-20
    partialflux=0;
end

end
