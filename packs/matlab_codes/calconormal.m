function [ai,aii]=calconormal(noi,noj,tensor,element)
global coord  centelem

% o conormal foi calculado utilizando a equação (3.1) e (3.2)
% do artigo Zhiming Gao and Jiming Wu 2015.
R=[0 1 0; -1 0 0;0 0 0];
% coordenada da face
IJ=coord(noj,:)-coord(noi,:);
% calculo da norma da face
normIJ=norm(IJ);
% calculo as coordenadas dos nos 
d=coord(noi,:); % primeiro vertice da face em questão
c=coord(noj,:); % segundo vertice da face em questão

%%
bi=dot(d-centelem(element,:),R*(c-centelem(element,:))');

ai= (R*(IJ/normIJ)')'*tensor*(R*(c-centelem(element,:))')/bi;
%%
bii=dot(c-centelem(element,:),R*(d-centelem(element,:))');

aii=(R*(IJ/normIJ)')'*tensor*(R*(d-centelem(element,:))')/bii;


end