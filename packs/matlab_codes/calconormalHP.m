function [alfai,alfaj]=calconormalHP(noi,noj,K,element,IJ,harmonicpoint)
global centelem
% do artigo Zhiming Gao and Jiming Wu 2015.

% alfai, alfaj --> coeficientes da combina��o linear, 3.1
% matriz de rota��o
R=[0 1 0; -1 0 0;0 0 0];
% coordenada da face
%IJ=coord(noj,:)-coord(noi,:);
% calculo da norma da face
normIJ=norm(IJ);

% calculo as coordenadas dos nos
d=harmonicpoint(noi,:); % primeiro vertice da face em quest�o
c=harmonicpoint(noj,:); % segundo vertice da face em quest�o

% baricentro do elemento em quest�o
center=centelem(element,:);
%%
bi=dot(d-center,R*(c-center)');

alfai= (R*(IJ/normIJ)')'*K*(R*(c-center)')/bi;
%%
bii=dot(c-center,R*(d-center)');

alfaj=(R*(IJ/normIJ)')'*K*(R*(d-center)')/bii;

end