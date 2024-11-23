import numpy as np
from sympy import *
from sympy import Symbol, Derivative, solve
from scipy.interpolate import InterpolatedUnivariateSpline


def inver(s):
    return (-(1250*s - 81*((10000*s)/81 + 4)**(1/2) + 162)/(1250*s))**(1/2)/2 + 1/2

def f(s):
    return 0.1296*(s**2/(s**2+(1-s)**2))

def df(s):
    return  0.1296*(-2*(s-1)*s/(2*s**2-2*s+1)**2);

def BL_exact(S0,x, fi, A, Swr, Sor, ftime):

    SL = 1 - Sor;
    SR = 0 + Swr;
    inf_point=0.5;
    St = 0.7071;

    tstep = 0.2;
    S = S0;
    t = ftime;
    for i in range(len(x)-1):
        if SR<SL and SL<=inf_point:
            vel = (A/fi)*(f(S[i])-f(S[i+1]))/(S[i]-S[i+1])
            pos = vel*t
            if x[i]<=pos:
                S[i]=S[i];
            elif x[i]>pos:
                S[i]=S[i+1];


        elif inf_point<=SR and SR<SL:

            if x[i]/t<=df(S[i]):
                S[i]=S[i]
            elif df(S[i])<x[i]/t and x[i]/t<df(S[i+1]):
                S[i]=inver(x[i]/t)
            elif df(S[i+1])<x[i]/t:
                S = S[i+1]

        elif SR<inf_point and inf_point<St and St<SL:
            ss = (f(St)-f(S[i+1]))/(St-S[i+1]);
            if (x[i]/t)<=df(S[i]):
                S[i]=S[i]
            elif df(S[i])<(x[i]/t) and ss>(x[i]/t):
                S[i]=inver(x[i]/t)
            elif ss<(x[i]/t):
                S[i]=S[i+1]

    return S
