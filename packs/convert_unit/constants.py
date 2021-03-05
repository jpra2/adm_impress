
k0 = 6894.76
k1 = 1e-3
k2 = 0.158987
k3 = 86400
k4 = 9.869233e-13
k5 = 3.28084
k6 = 12
k7 = 9.869233e-16
k8 = 0.158987/(86400)

def psi_to_Pa(val=1.0):
    '''
    psi em pascal
    '''
    return k0*val

def Pa_to_psi(val=1.0):
    """
    pascal em psi
    """
    return val/k0

def cp_to_Pas(val=1.0):
    '''
    centipoise em pascal segundo
    '''
    return val*k1

def Pas_to_cp(val):
    '''
    pascal segundo em centipoise
    '''
    return val/k1

def bbl_to_m3(val=1.0):
    '''
    bbl em metro cubico
    '''
    return val*k2

def m3_to_bbl(val=1.0):
    '''
    metro cubico em bbl
    '''
    return val/k2

def d_to_s(val=1.0):
    """
    dia em segundo
    """
    return val*k3

def s_to_d(val=1.0):
    """
    segundo em dia
    """
    return val/k3

def darcy_to_m2(val=1.0):
    '''
    darcy to m2
    '''
    return val*k4

def m2_to_darcy(val=1.0):
    '''
     m2 to darcy
    '''
    return val/k4

def m_to_pe(val=1.0):
    '''
    metro para pe
    '''
    return val*k5

def pe_to_m(val=1.0):
    '''
    pe to metro
    '''
    return val/k5

def pe_to_pol(val=1.0):
    '''
    pe to polegada
    '''
    return val*k6

def pol_to_pe(val=1.0):
    '''
    polegada em pe
    '''
    return val/k6

def milidarcy_to_m2(val=1.0):
    '''
    milidarcy em metro quadrado
    '''
    return val*k7

def m2_to_milidarcy(val=1.0):
    '''
    metro quadrado em milidarcy
    '''
    return val/k7

def bbldia_to_m3seg(val=1.0):
    '''
    bbl/dia to m3/s
    '''
    return val*k8

def m3seg_to_bbldia(val=1.0):
    '''
    m3/s to bbl/dia
    '''
    return val/k8
