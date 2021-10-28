import numpy as np
from sympy import *
from sympy import Symbol, Derivative, solve
from scipy.interpolate import InterpolatedUnivariateSpline

''' ------------------------------- Geometria ------------------------------ '''
''' Comprimento (m):'''
L = 300; #300

''' Area da secao transversal (m?):'''
A = 1;

''' ------------------ Propriedades dos fluidos e da rocha -----------------'''
''' Porosidade (%):'''
phi = 0.2;

''' Viscosidades (cp):'''
mio = 1e-3
miw = 1e-3

''' ------------------- Condi??es iniciais e de contorno -------------------'''
''' Satura??o de ?leo residual (%):(condicao de contorno ) '''
Sor = 0.;

''' Satura??o de ?gua (%) (condi??o inicial Srw): '''
Swr = 0.;

''' Vaz?o de inje??o de ?gua (m? std/d) (Como ? std tem que usar o Bw): '''
Q = 3.e-7*60*24*60;

''' Fator volume da forma??o da ?gua (m?/m? std):'''
Bw = 1

time = 1500 #d

''' Criando vari?veis simb?licas (para poder mudar o modelos): '''

'''Model parameters durlofsky'''
row = 1.000;#density of water
roo = 0.781855078;#density of oleo
g = 9.8;#gravity
''' Permeabilidades relativas (Usando o modelo de Corey): '''
no = 2;
nw = 2;

sw  = symbols('sw')

''' Satura??o normalizada: '''
swn = (sw-Swr)/(1-Swr-Sor); #ok<NODEF>
kro = (1-swn)**no;
krw = swn**nw;

''' Fluxo fracion?rio (lambda_w/lmabda_t): '''

fw = (krw/miw)/((krw/miw)+(kro/mio));

''' derivada simbolica da fun??o de fluxo: '''

dfwdsw = fw.diff(sw);
#p1 = InterpolatedUnDarlanariateSpline(sw,fw)
#p2 = InterpolatedUnivariateSpline.derivative(p1)
#dfwdsw = p2(sw)

''' Welge tangent saturation (Le Veque 2002) - Resolve a equa??o para achar
 a satura??o da frente:'''
sta = np.array(solve((sw-Swr)*dfwdsw-fw,sw))

st = sta[np.isreal(sta)];
roots = np.double(st);
roots = roots[roots<=1];
roots = roots[roots>=0];
swf = max(roots);

''' Satura??o da frente:'''
#sw = swf; #ok<NASGU>

''' Fun??o satura??o em fun??o da posi??o:
%
%           / Sw(i) = Sw(x) = (inversa de x(sw)), se x(i) <= xf
%   Sw(i) = |
%           \ Sw(i) = Swr, se x(i) >= xf

% Como ja podemos calcular o valor da posi??o que em determinado tempo ter?
% uma determinada satua??o criamos um vetor com espa?amento constante de
% satura??es, calculamos o x correspondente e plotamos (x,sw).'''




''' Considerando-se que a vaz?o ? constante no tempo e igual a Q a inje??o
% acumulada ? dada por: W = Q*t. Desta forma pode-se escrever a posi??o
% como um perfil (multiplo da derivada do fluxo) vezes o tempo
% Perfil:
% Equa??o para a posi??o (dividida pelo tempo), Eq. (14.214) do livro
% Engenharia de reservat?rios de petr?leo (Livro verde do Adalberto):'''

f = lambdify(sw,dfwdsw)
sw = np.linspace(Swr, 1-Sor, 1/0.001+1);

profile = (Bw*Q/(A*phi*L))*np.double(f(sw));
xp = profile*time;

xf = np.max(xp[sw>=swf])
Sw = np.ones_like(sw)*(-1)
Sw[sw>=swf] = sw[sw>=swf]
xp_xf = max(xp[xp==xf])
var = np.argwhere((xp==xp_xf))
Sw[max(var)] = swf
Sw[max(var)-1] = Swr
Sw[max(var)-2] = Swr

xD = np.ones_like(xp)*(-1)
xD[sw>=swf] = xp[sw>=swf]
xD[max(var)] = xf
xD[max(var)-1] = xf
xD[max(var)-2] = 1
x_ans = xD[max(var)[0]-2:]
Sw_ans = Sw[max(var)[0]-2:]

import pdb; pdb.set_trace()
with open('Sw_BL_FR_semi_analytical.txt', 'w') as f:
    for item in Sw_ans:
        f.write("%s\n" % item)

with open('x_BL_FR_semi_analytical.txt', 'w') as f:
    for item in x_ans:
        f.write("%s\n" % item)
import pdb; pdb.set_trace()
