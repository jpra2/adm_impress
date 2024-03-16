function [elemphant,elembedge,face_bedge,peso,nu]=premultidimensional(N)

%% utilize somente para Método Multidimensional 
peso=0.1;
nu=0.75;

%% elementos no contorno no método Multidimensional
[elemphant,elembedge,face_bedge]=elemphanton(N);

end