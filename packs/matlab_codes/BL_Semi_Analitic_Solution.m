%Equa�ao de Buckley-Leverett

%   15-09-2006
%   Joao Paulo Silva Queiroz

%   Este programa constroi o perfil de saturacao de um escoamento bifasico,
%   uinidimensional (oleo e agua), utilizando a equacao de Buckley-Leverett e
%   o metodo de Welge para sua resolucao.


% clc;

%   Dados de entrada
npoin=170;
Swr=0.9;     %Satura��o residual da �gua
Sor=0.1;     %Satura��o residual do �leo
mi_w=1;      %Viscosiade da �gua
mi_o=1;      %Viscosiade do �leo
A=1;         %Area normal ao escoamento
%fi=.2;       %Porosidade do meio
fi=1;
q=1;
%q=.2;        %Taxa de inje�ao de �gua
t=0.2;        %Tempo transcorrido
L=1;         %Dist�ncia entre injetor e produtor


%   Satura�ao no choque
%   M�todo das secantes

%   Valores iniciais para o metodo
Swf=1-Sor;
Sw0=Swf/2;

while abs(Swf-Sw0)>1e-12
    Sw=Swf;
    
%   Calculo dos pontos da funcao de fluxo fracional
    fw0=1/(1+mi_o/mi_w*((1-Sw0-Swr)/(Sw0-Swr))^2);
    fw=1/(1+mi_o/mi_w*((1-Sw-Swr)/(Sw-Swr))^2);
    
%   Calculo da derivada da funcao de fluxo fracional
    dfw_dSw0=-2*mi_w/mi_o*((1-Sw0-Swr)/(Sw0-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Sw0-Swr)/(Sw0-Swr))^2)^2;
    dfw_dSw=-2*mi_w/mi_o*((1-Sw-Swr)/(Sw-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Sw-Swr)/(Sw-Swr))^2)^2;

%   Equa�ao da reta tangente    
    y0=dfw_dSw0*(Sw0-Swr);
    y=dfw_dSw*(Sw-Swr);
    
%   Swf=> Zero da fun�ao g=y-fw
    g0=y0-fw0;
    g=y-fw;
    Swf=Sw-(Sw-Sw0)*g/(g-g0);
    Sw0=Sw;
    
end;

%   Eixo x

delt_x=L/(npoin-1);
x=0:delt_x:L;

%   Posicao do choque
dfw_dSw_Swf=-2*mi_w/mi_o*((1-Swf-Swr)/(Swf-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Swf-Swr)/(Swf-Swr))^2)^2;
xf=q*t/(A*fi)*dfw_dSw_Swf;

Sw(2:length(x))=0.5;
Sw(1)=1-Sor;
ii=1;
while x(ii)<=xf
    ii=ii+1;
    Sw0=.9;
    while abs(Sw(ii)-Sw0)>1e-12
%         disp(ii);
%       Calculo da derivada da funcao de fluxo fracional
        dfw_dSw0=-2*mi_w/mi_o*((1-Sw0-Swr)/(Sw0-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Sw0-Swr)/(Sw0-Swr))^2)^2;
        dfw_dSw=-2*mi_w/mi_o*((1-Sw(ii)-Swr)/(Sw(ii)-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Sw(ii)-Swr)/(Sw(ii)-Swr))^2)^2;
%       Saturacao em cada x
        Sw1=Sw(ii)-(Sw(ii)-Sw0)*(dfw_dSw-A*fi/(q*t)*x(ii))/(dfw_dSw-dfw_dSw0);
        Sw0=Sw(ii);
        Sw(ii)=Sw1;
    end;
end;
Sw(ii:length(x))=Swr;


%   Grafico

hold all;
plot(x,Sw,'r-');
