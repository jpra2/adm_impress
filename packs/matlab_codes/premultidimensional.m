function [elemphant,elembedge,face_bedge,peso,nu]=premultidimensional(N)

%% utilize somente para M�todo Multidimensional 
peso=0.1;
nu=0.75;

%% elementos no contorno no m�todo Multidimensional
[elemphant,elembedge,face_bedge]=elemphanton(N);

end