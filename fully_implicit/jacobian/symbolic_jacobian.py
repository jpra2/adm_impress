import sympy as sym
import yaml

class symbolic_J:
    def __init__(self):
        with open('input_cards/fluid_properties.yml', 'r') as f:
            fluid_properties = yaml.safe_load(f)

        Swc=fluid_properties['Swc'] #Saturação residual da água
        Sor=fluid_properties['Sor'] #Saturação residual do óleo
        n_o=fluid_properties['n_o'] #Coeficiente de Corey brooks do óleo
        n_w=fluid_properties['n_w'] #Coeficiente de Corey brooks da água
        r_o=fluid_properties['gama_o'] #Densidade do óleo
        r_w=fluid_properties['gama_w'] #Densidade da água
        mi_o=fluid_properties['mi_o'] #Viscosidade cinemática da água
        mi_w=fluid_properties['mi_w'] #Viscosidade cinemática do óleo
        relative_permeability_model=fluid_properties["relative_permeability"]

        T, S_up, Sw, So, Swn, Son, Dt, k, phi, p_i, p_j, Dx, Dy=sym.symbols("T S Sw So Swn Son Dt k phi p_i p_j Dx Dy")
        if relative_permeability_model=='BrooksAndCorey':
            krw=((Sw - Swc)/(1 - Swc - Sor))**n_o
            kro=(1-(Sw - Swc)/(1 - Swc - Sor))**n_w

        lam_w=krw/mi_w
        lam_o=kro/mi_o

        self.F_w=r_w*lam_w*T*(p_j-p_i)
        self.F_o=r_o*lam_o*T*(p_j-p_i)

        self.acum_w=Dx*Dy*phi*r_w*(Sw-Swn)/Dt
        self.acum_o=Dx*Dy*phi*r_o*(1-Sw-(1-Swn))/Dt

        J=[[sym.diff(self.F_o,p_i), sym.diff(self.F_o,Sw)],[sym.diff(self.F_w,p_i), sym.diff(self.F_w,Sw)]]

        self.J=J
        self.c_o=sym.diff(self.acum_o,Sw)
        self.c_w=sym.diff(self.acum_w,Sw)
        # import pdb; pdb.set_trace()
