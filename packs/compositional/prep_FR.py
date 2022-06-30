import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from sympy import Symbol, diff, poly
import quadpy
from scipy.special import legendre, jacobi, gamma


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

    #dgLB_lump_Lo =  polynomial_der_coeffs(g2_coefs_LB.c[::-1])[:,:-1]
    #dgRB_lump_Lo =  polynomial_der_coeffs(g2_coefs_RB.c[::-1])[:,:-1]
    # starts with lower power term, and ends with a minor order than its original
    # function, i.e. dgLB[:,-1] is not zero because the vector has the number of
    # elements equal to the derivative order. This makes easier when computing dFk
    return dgRB_Radau, dgLB_Radau, LRadau_coefs.c[::-1], RRadau_coefs.c[::-1]
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

def dJacobiP(order, alpha, beta, points):
    return  np.sqrt(order*(order+alpha+beta+1))*jacobiP_((order-1)*np.sign(order), \
        alpha+1, beta+1, points)

def jacobiP_(order, alpha, beta, points):
    xp = points; dims = np.size(points);
    #if (dims(2)==1); xp = xp'; end;
    PL = np.zeros((order+1,np.size(points)));

    #Initial values P_0(x) and P_1(x)

    gamma0 = 2.**(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)* \
            gamma(beta+1)/gamma(alpha+beta+1)
    PL[0,:] = 1.0/np.sqrt(gamma0)
    if (order==0): return PL.T

    gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
    PL[1,:] = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/np.sqrt(gamma1)
    if (order==1): return PL[order,:].T

    #Repeat value in recurrence.
    aold = 2/(2+alpha+beta)*np.sqrt((alpha+1)*(beta+1)/(alpha+beta+3))

    #Forward recurrence using the symmetry of the recurrence.
    for i in range(1,order):
        h1 = 2*(i)+alpha+beta
        anew = 2/(h1+2)*np.sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)* \
            (i+1+beta)/(h1+1)/(h1+3))
        bnew = - (alpha**2-beta**2)/h1/(h1+2)
        PL[i+1,:] = 1/anew*( -aold*PL[i-1,:] + (xp-bnew)*PL[i,:])
        aold =anew

    return PL[order,:].T

def jacobiP(order, alpha, beta, points):
    return jacobi(order, alpha, beta) * np.sqrt((2*order+1)/2)

def gradVandermonde(n_points, points):
    Vr = np.zeros((n_points,n_points)) # Allocate

    for i in range(n_points):	# All Polynomial Degrees up to kDeg
        for j in range(n_points):
            Vr[j,i] = dJacobiP(i,0,0,points[j])
    return Vr

def Vandermonde2(n_points, points):
    Vr = np.zeros((n_points,n_points)) # Allocate
    for i in range(n_points):	# All Polynomial Degrees up to kDeg
        for j in range(n_points):
            Vr[j,i] = jacobiP_(i,0,0,points[j])
    return Vr

def run(M):
    global n_points
    global points
    global weights
    global L
    global dL
    global dgRB
    global dgLB
    global V
    global x_points
    global v0
    global vols_vec
    global Dr

    n_points = data_loaded['compositional_data']['FR']['order']
    GL = quadpy.c1.gauss_lobatto(n_points)
    points = GL.points
    weights = GL.weights
    dgRB, dgLB, gRB, gLB = correction_function(n_points)
    L, dL = Lagrange_poly(n_points, points)

    points[abs(points)<1e-15] = 0
    L[abs(L)<5e-15] = 0
    dL[abs(dL)<5e-15] = 0
    '''points = np.round(points,2)
    L = np.round(L,2)
    dL = np.round(dL,2)'''


    V = np.polynomial.legendre.legvander(points,n_points-1)
    x_points, v0, vols_vec = auxiliary_terms(M, points, n_points)

    #V_H = Vandermonde2(n_points, points)
    #GV_H = gradVandermonde(n_points, points)
    #Dr = GV_H@np.linalg.inv(V_H)
