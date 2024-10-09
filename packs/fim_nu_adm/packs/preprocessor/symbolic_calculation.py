from ..inputs import finescale_inputs as inputs
import sympy as sym
import yaml

class symbolic_J:
    def __init__(self):
        fluid_properties = inputs['fluid_properties']
        Swc=fluid_properties['Swc'] #Saturação residual da água
        Sor=fluid_properties['Sor'] #Saturação residual do óleo
        n_o=fluid_properties['n_o'] #Coeficiente de Corey brooks do óleo
        n_w=fluid_properties['n_w'] #Coeficiente de Corey brooks da água

        r_o=fluid_properties['gama_o'] #Densidade do óleo
        r_w=fluid_properties['gama_w'] #Densidade da água
        mi_o=fluid_properties['mi_o'] #Viscosidade cinemática da água
        mi_w=fluid_properties['mi_w'] #Viscosidade cinemática do óleo
        relative_permeability_model=fluid_properties["relative_permeability"]

        T, S_up, Sw, So, Swn, Son, Dt, k, p_i, p_j, Vp, gh=sym.symbols("T S Sw So Swn Son Dt k p_i p_j Vp gh")
        if relative_permeability_model=='BrooksAndCorey':
            krw=((Sw - Swc)/(1 - Swc - Sor))**n_o
            kro=(1-(Sw - Swc)/(1 - Swc - Sor))**n_w

        lam_w=krw/mi_w
        lam_o=kro/mi_o
        self.g=fluid_properties['g']
        self.r_o=r_o
        self.r_w=r_w
        self.F_w=r_w*lam_w*T*(p_j-p_i+r_w*gh)
        self.F_o=r_o*lam_o*T*(p_j-p_i+r_o*gh)

        self.acum_w=Vp*r_w*(Sw-Swn)/Dt
        self.acum_o=Vp*r_o*(1-Sw-(1-Swn))/Dt

        J=[[sym.diff(self.F_o,p_i), sym.diff(self.F_o,Sw)],[sym.diff(self.F_w,p_i), sym.diff(self.F_w,Sw)]]

        self.J=J

        self.c_o=sym.diff(self.acum_o,Sw)
        self.c_w=sym.diff(self.acum_w,Sw)
        
        # import pdb; pdb.set_trace()

        self.c_o=(sym.lambdify((Vp,Dt),self.c_o))
        self.c_w=(sym.lambdify((Vp,Dt),self.c_w))
        self.acum_o = sym.lambdify((Vp,Dt,Sw,Swn),self.acum_o)
        self.acum_w = sym.lambdify((Vp,Dt,Sw,Swn),self.acum_w)
        # import pdb; pdb.set_trace()
        self.J[0][0] = sym.lambdify((T,Sw),self.J[0][0])
        self.J[0][1] = sym.lambdify((T,Sw, p_i, p_j, gh),self.J[0][1])
        self.J[1][0] = sym.lambdify((T,Sw),self.J[1][0])
        self.J[1][1] = sym.lambdify((T,Sw, p_i, p_j, gh),self.J[1][1])
        self.F_o = sym.lambdify((T,Sw, p_i, p_j, gh),self.F_o)
        self.F_w = sym.lambdify((T,Sw, p_i, p_j, gh),self.F_w)
