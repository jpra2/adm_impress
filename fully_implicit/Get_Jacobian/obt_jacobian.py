import sympy as sym
# from sympy import Heaviside as H

class F_Jacobian:
    def __init__(self):
        self.teste=5.0
        Swc=0.2 #Saturação residual da água
        Sor=0.2 #Saturação residual do óleo
        n_o=2.0 #Coeficiente de Corey brooks do óleo
        n_w=2.0 #Coeficiente de Corey brooks da água
        r_o=0.8 #Densidade do óleo
        r_w=1.0 #Densidade da água
        mi_o=5.0 #Viscosidade cinemática da água
        mi_w=1.0 #Viscosidade cinemática do óleo

        T, S_up, Sw, So, Swn, Son, Dt, k, phi, p_i, p_j, Dx, Dy=sym.symbols("T S Sw So Swn Son Dt k phi p_i p_j Dx Dy")

        krw=((Sw - Swc)/(1 - Swc - Sor))**n_o
        kro=(1-(Sw - Swc)/(1 - Swc - Sor))**n_w

        lam_w=k*krw/mi_w
        lam_o=k*kro/mi_o

        self.F_w=r_w*lam_w*T*(p_j-p_i)
        self.F_o=r_o*lam_o*T*(p_j-p_i)

        self.acum_w=Dx*Dy*phi*r_w*(Sw-Swn)/Dt
        self.acum_o=Dx*Dy*phi*r_o*(1-Sw-(1-Swn))/Dt

        J=[[sym.diff(self.F_o,p_i), sym.diff(self.F_o,Sw)],[sym.diff(self.F_w,p_i), sym.diff(self.F_w,Sw)]]

        self.J=J
        self.c_o=sym.diff(self.acum_o,Sw)
        self.c_w=sym.diff(self.acum_w,Sw)
