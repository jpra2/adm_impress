import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from sympy import Symbol, diff, poly
import quadpy
from scipy.special import legendre


global x
x = Symbol('x')

def Lagrange_poly(n_points, points):
    np.seterr(divide='ignore', invalid='ignore')
    L = np.empty((n_points, n_points))
    dL = np.empty_like(L)
    for i in range(n_points):
        L_poly = 1
        for j in range(n_points):
            if j!=i:
                L_poly = L_poly * (x - points[j]) * \
                (1 / (points[i] - points[j]))
        L[i,:] = poly(L_poly,x).coeffs() #starting from higher power term
        L[i,:] = L[i,::-1] #reversing order
        dL[i,:] = polynomial_der_coeffs(L[i])
    return L, dL

def correction_function(n_points):
    'Radau coefficients'
    g_coefs_LegendreK = legendre(n_points)
    g_coefs_LegendreK_1 = np.zeros_like(g_coefs_LegendreK)
    g_coefs_LegendreK_1[1:] = legendre(n_points-1) #np.polynomial.legendre.legfromroots(pointsK_1)
    RRadau_coefs = (-1) ** (n_points)/2 * (g_coefs_LegendreK - g_coefs_LegendreK_1)
    LRadau_coefs = 1/2 * (g_coefs_LegendreK + g_coefs_LegendreK_1)

    dgLB_Radau =  polynomial_der_coeffs(RRadau_coefs.c[::-1])[:,:-1]
    dgRB_Radau =  polynomial_der_coeffs(LRadau_coefs.c[::-1])[:,:-1]

    'Lump correction_function'
    g_coefs_LegendreK_2 = np.zeros_like(g_coefs_LegendreK)
    g_coefs_LegendreK_2[2:] = legendre(n_points-2)
    RRadau_coefs_K_1 = (-1) ** (n_points-1)/2 * (g_coefs_LegendreK_1 - g_coefs_LegendreK_2)
    LRadau_coefs_K_1 =  1/2 * (g_coefs_LegendreK_1 + g_coefs_LegendreK_2)
    g2_coefs_LB = (n_points-1)/(2*n_points-1) * RRadau_coefs + n_points/(2*n_points-1) * RRadau_coefs_K_1
    g2_coefs_RB = (n_points-1)/(2*n_points-1) * LRadau_coefs + n_points/(2*n_points-1) * LRadau_coefs_K_1

    dgLB_lump_Lo =  polynomial_der_coeffs(g2_coefs_LB.c[::-1])[:,:-1]
    dgRB_lump_Lo =  polynomial_der_coeffs(g2_coefs_RB.c[::-1])[:,:-1]
    # starts with lower power term, and ends with a minor order than its original
    # function, i.e. dgLB[:,-1] is not zero because the vector has the number of
    # elements equal to the derivative order. This makes easier when computing dFk
    return dgRB_Radau, dgLB_Radau, RRadau_coefs.c[::-1], LRadau_coefs.c[::-1]
    #return dgRB_lump_Lo, dgLB_lump_Lo, g2_coefs_RB.c[::-1], g2_coefs_LB.c[::-1]

def polynomial_der_coeffs(P):
    'get the derivative from a given polynomial P given in coefficient form \
    where the first column corresponds to the lowest order polynomial coeff. \
    This function returns the derivative coefficients with the number of columns \
    corresponding to the order of the vector P, to be derivated. P can also be a matrix'
    dP = np.array([[c*k for k,c in enumerate(P)]])
    dP[:,0:-1] = dP[:,1:]
    dP[:,-1] = 0
    return dP

def RT0_shape_functions():
    phi = np.empty((n_points,2))
    phi[:,0] = 1 / 4 * (1 + x)
    phi[:,1] = 1 / 4 * (1 - x)
    return phi

def auxiliary_terms(M, points, n_points):
    x_points = np.array([points**i for i in range(n_points)])
    c_int = M.faces.center(M.faces.internal)
    c_vols = M.volumes.center(M.volumes.all)
    pos = (c_int[:,np.newaxis,:] - c_vols[ctes.v0]).sum(axis=2)
    v0 = np.copy(ctes.v0)
    v0[:,0] = ctes.v0[pos>0]
    v0[:,1] = ctes.v0[pos<0]

    'Get neigboring cells values'
    vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
    lines = np.arange(ctes.n_internal_faces)
    vols_vec[v0[:,0],1] = lines
    vols_vec[v0[:,1],0] = lines
    contour_vols = np.argwhere(vols_vec<0)[:,0]
    vols_vec[vols_vec < 0] = vols_vec[contour_vols,:][vols_vec[contour_vols]>=0]

    return x_points, v0, vols_vec

def run(M):
    global n_points
    global points
    global weights
    global L
    global dL
    global dgRB, gRB
    global dgLB, gLB
    global V
    global x_points
    global v0
    global vols_vec

    n_points = data_loaded['compositional_data']['FR']['order']
    GL = quadpy.c1.gauss_lobatto(n_points)
    points = GL.points
    weights = GL.weights
    dgRB, dgLB, gRB, gLB = correction_function(n_points)
    L, dL = Lagrange_poly(n_points, points)
    #points = np.round(points,2)
    #L = np.round(L,2)
    #dL = np.round(dL,2)
    V = np.polynomial.legendre.legvander(points,n_points-1)
    x_points, v0, vols_vec = auxiliary_terms(M, points, n_points)
