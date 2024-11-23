
k0 = 6894.76
k1 = 1e-3
k2 = 0.158987
k3 = 86400
k4 = 9.869233e-13
k5 = 3.28084
k6 = 12
k7 = 9.869233e-16
k8 = 0.158987/(86400)

import pint
import os
path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, 'const_definitions.txt')
ureg = pint.UnitRegistry()
ureg.load_definitions(path)

def psi_to_Pa(val=1.0):
    '''
    psi em pascal
    '''

    psi = val*ureg.psi
    pa = psi.to(ureg.Pa).magnitude

    return pa

def Pa_to_psi(val=1.0):
    """
    pascal em psi
    """
    pa = val*ureg.Pa
    psi = pa.to(ureg.psi).magnitude

    return psi

def cp_to_Pas(val=1.0):
    '''
    centipoise em pascal segundo
    '''
    cp = val*ureg.centipoise
    Pas = cp.to(ureg.Pa * ureg.s).magnitude

    return Pas

def Pas_to_cp(val):
    '''
    pascal segundo em centipoise
    '''
    Pas = val*(ureg.Pa * ureg.s)
    cp = Pas.to(ureg.centipoise).magnitude
    return cp

def bbl_to_m3(val=1.0):
    '''
    bbl em metro cubico
    '''
    bbl = val*ureg.bbl
    m3 = bbl.to(ureg.m**3).magnitude
    return m3

def m3_to_bbl(val=1.0):
    '''
    metro cubico em bbl
    '''
    m3 = val*(ureg.m**3)
    bbl = m3.to(ureg.bbl).magnitude
    return bbl

def d_to_s(val=1.0):
    """
    dia em segundo
    """
    d = val*ureg.day
    seg = d.to(ureg.s).magnitude
    return seg

def s_to_d(val=1.0):
    """
    segundo em dia
    """
    seg = val*ureg.s
    d = seg.to(ureg.day).magnitude
    return d

def darcy_to_m2(val=1.0):
    '''
    darcy to m2
    '''
    darcy = val*ureg.darcy
    m2 = darcy.to(ureg.m**2).magnitude
    return m2

def m2_to_darcy(val=1.0):
    '''
     m2 to darcy
    '''
    m2 = val*(ureg.m**2)
    darcy = m2.to(ureg.darcy).magnitude
    return darcy
    # return val/k4

def m_to_pe(val=1.0):
    '''
    metro para pe
    '''
    m = val*ureg.m
    ft = m.to(ureg.ft).magnitude
    return ft
    # return val*k5

def pe_to_m(val=1.0):
    '''
    pe to metro
    '''
    ft = val*ureg.ft
    m = ft.to(ureg.m).magnitude
    return m
    # return val/k5

def pe_to_pol(val=1.0):
    '''
    pe to polegada
    '''
    ft = val*ureg.ft
    pol = ft.to(ureg.inch).magnitude
    return pol
    # return val*k6

def pol_to_pe(val=1.0):
    '''
    polegada em pe
    '''
    pol = val*ureg.inch
    ft = pol.to(ureg.ft).magnitude
    return ft
    # return val/k6

def milidarcy_to_m2(val=1.0):
    '''
    milidarcy em metro quadrado
    '''
    md = val*ureg.md
    m2 = md.to(ureg.m**2).magnitude
    return m2
    # return val*k7

def m2_to_milidarcy(val=1.0):
    '''
    metro quadrado em milidarcy
    '''
    m2 = val*(ureg.m**2)
    md = m2.to(ureg.md).magnitude
    return md
    # return val/k7

def bbldia_to_m3seg(val=1.0):
    '''
    bbl/dia to m3/s
    '''
    bbl_dia = val*(ureg.bbl / ureg.day)
    m3s = bbl_dia.to((ureg.m**3) / ureg.s).magnitude
    return m3s
    # return val*k8

def m3seg_to_bbldia(val=1.0):
    '''
    m3/s to bbl/dia
    '''
    m3s = val*((ureg.m**3) / ureg.s)
    bbldia = m3s.to(ureg.bbl / ureg.day).magnitude
    return bbldia
    # return val/k8

def atm_to_psi(val=1.0):

    atm = val*ureg.atm
    psi = atm.to(ureg.psi).magnitude
    return psi

def psi_to_atm(val=1.0):
    psi = val*ureg.psi
    atm = psi.to(ureg.atm).magnitude
    return atm

def atm_to_pa(val=1.0):
    atm = val*(ureg.atm)
    pa = atm.to(ureg.Pa).magnitude
    return pa

def pa_to_atm(val=1.0):
    Pa = val*(ureg.Pa)
    atm = Pa.to(ureg.atm).magnitude
    return atm
