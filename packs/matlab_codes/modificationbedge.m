% algumas linhas de esta funcao deve ser ativado para malhas bem
% especificas
function [bedge]=modificationbedge(bedge)

%%
% bedge(217,4)=102;
% bedge(223,4)=102;
% bedge(229,4)=102;
% bedge(235,4)=102;
%% só em malha com furo 1
% x=bedge(217:240,1);
% y=bedge(217:240,2);
% bedge(217:240,1)=y;
% bedge(217:240,2)=x;
%% só em malha com furo 2
% x=bedge(129:144,1);
% y=bedge(129:144,2);
% bedge(129:144,1)=y;
% bedge(129:144,2)=x;
%% distorcao de malhas estruturadas
% esta funcao pode ser ativado se deseja distorcer alguma malha estruturada
%[auxcoord]=distortedramd;

%% somente use nas malhas "Tipo1malha1", "Tipo1malha2", "Tipo1malha3" e "Tipo1malha4"
%Tipo4Malha1_1, Tipo4Malha1_2, Tipo4Malha1_3, Tipo4Malha1_4
% Historial para malha "Tipo1malha0" não habilite nada, el ya viene
% ordenado en sentido antihorario tanto en elemento como no contorno
% para malha "Tipo1malha1" "Tipo1malha2 " "Tipo1malha3" e "Tipo1malha4"
% vamos habilitar na linha 803-806 do preprocessador,
% e na linha, caso contrario vai dar erro no calculo do esurn1 esurn2 e no
% nsurn1 e no nsurn2

%   x=bedge(:,1);
%   y=bedge(:,2);
%   bedge(:,1)=y;
%   bedge(:,2)=x;
%   x1=elem(:,1);
%   x2=elem(:,3);
%   elem(:,1)=x2;
%   elem(:,3)=x1;

%% Modificação Malha Kershaw
%bedge(:,4:5)=101;
%----------------------------
%% tratamento malha Hermeline
%bedge(:,4:5)=101;
% malha 16x16
%x=bedge(16:24,1);
%y=bedge(16:24,2);
%bedge(16:24,1)=y;
%bedge(16:24,2)=x;
% malha 32x32
% x=bedge(33:48,1);
% y=bedge(33:48,2);
% bedge(33:48,1)=y;
% bedge(33:48,2)=x;
% malha 64x64
% x=bedge(64:96,1);
% y=bedge(64:96,2);
% bedge(64:96,1)=y;
% bedge(64:96,2)=x;
% malha 128x128
%x=bedge(128:192,1);
%y=bedge(128:192,2);
%bedge(128:192,1)=y;
%bedge(128:192,2)=x;

end