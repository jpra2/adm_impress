import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from sympy import Symbol, diff, poly
#import quadpy
#from quadpy.c1 import gauss_lobatto
from scipy.special import legendre, jacobi, gamma
from scipy.interpolate import lagrange

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

        if len(poly(L_poly,x).coeffs())<n_points:
            L[i,[-1,1]] = poly(L_poly,x).coeffs()
        else: L[i,:] = poly(L_poly,x).coeffs() #starting from higher power term
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

def RT0_shape_functions(n_points, points):
    phi = np.empty((n_points,2))
    phi[:,0] = 1 / 4 * (1 + points)
    phi[:,1] = 1 / 4 * (1 - points)
    return phi

def auxiliary_terms_2D(faces_around_vols_2D_idx):
    contour_faces_around_vols_2D = ~np.isin(faces_around_vols_2D_idx,ctes.internal_faces)
    vols_vec = np.copy(faces_around_vols_2D_idx)
    lines = np.arange(ctes.n_internal_faces)
    int_faces = ctes.internal_faces
    vols_vec[contour_faces_around_vols_2D] = -1
    rr=vols_vec[vols_vec>=0]
    vols_vec[vols_vec>=0] = np.searchsorted(ctes.internal_faces,rr)
    bool = vols_vec<0
    pares = {(0,2),(1,3),(2,0),(3,1)}
    for par in pares:
        vols_vec[bool[:,par[0]],par[0]] = vols_vec[bool[:,par[0]],par[1]]
    return vols_vec

def auxiliary_terms(faces_around_vols_2D_idx):
    x_points = np.array([points**i for i in range(n_points)])
    lines = np.arange(ctes.n_internal_faces)
    'Get neigboring cells values -  this only works for 1D problems'
    vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
    lines = np.arange(ctes.n_internal_faces)
    vols_vec[ctes.v0[:,0],1] = lines
    vols_vec[ctes.v0[:,1],0] = lines
    contour_vols = np.argwhere(vols_vec<0)[:,0]
    vols_vec[vols_vec < 0] = vols_vec[contour_vols,:][vols_vec[contour_vols]>=0]
    return x_points, vols_vec

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

def bilinear_shape_functions_2D(n_points,x_points, y_points):
    #function to map from the reference domain to the physical one
    Ns = np.empty((n_points,4))
    Ns[:,0] = 1/4*(1-x_points)*(1-y_points)
    Ns[:,1] = 1/4*(1+x_points)*(1-y_points)
    Ns[:,2] = 1/4*(1+x_points)*(1+y_points)
    Ns[:,3] = 1/4*(1-x_points)*(1+y_points)
    return Ns

def map_faces_and_points_in_vols(M, Ns):
    faces_vols = M.volumes.bridge_adjacencies(M.volumes.all,3,2)
    faces_vols_center_fl = M.faces.center(faces_vols.flatten())
    faces_vols_center = np.concatenate(np.split(faces_vols_center_fl[:,np.newaxis,:],len(faces_vols_center_fl)/6),axis=1)
    ff=faces_vols_center.transpose(1,0,2)
    faces_vols_2D_indx = faces_vols[ff[:,:,-1]==0] #faces index that represent the volumes in 2D
    nodes_faces_vols_2D =  M.data['faces_nodes'][faces_vols_2D_indx]
    nodes_faces_vols_2D_coords_fl = M.nodes.center(nodes_faces_vols_2D.flatten())
    nodes_faces_vols_2D_coords = np.concatenate(np.split(nodes_faces_vols_2D_coords_fl[:,np.newaxis,:],len(nodes_faces_vols_2D_coords_fl[:,0])/4),axis=1)
    nodes_faces_vols_2D_coords = nodes_faces_vols_2D_coords.transpose(1,0,2)[...,:2]

    v0 = nodes_faces_vols_2D_coords[:,0]

    angles = np.empty((ctes.n_volumes, 3))

    for i in range(1,4):
        v = nodes_faces_vols_2D_coords[:,i] - v0
        angle = np.arctan(v[:,1]/v[:,0])
        angle[(v[:,0]<0)*(angle==0)] = 2*np.pi # correcting for negative x value (it was 0 degrees when it should be 180)
        angle[(v[:,0]<0)] -= np.pi # correcting for negative x value
        angle[angle<0]+=2*np.pi*np.ones_like(angle[angle<0]) # map the angle to 0,2pi interval
        angles[:,i-1] = angle
        if i==2:# to correct the loop direction (right hand rule)
            angles[angles[:,i-1]>=(angles[:,i-2]+np.pi),0] += 2*np.pi
            angles[angles[:,i-2]>=(angles[:,i-1]+np.pi),1] += 2*np.pi
    angles_argsort = np.argsort(angles)

    vs = nodes_faces_vols_2D_coords[:,1:]
    vs_reord_aux = vs[:,angles_argsort] #rearrange vectors in the correct order,
    # however this gets too much info just to vectorize.
    #The info we need is in the diagonal:
    vs_reord = np.diagonal(vs_reord_aux)
    vs_reord = vs_reord.transpose(2,0,1) #correct matrix shape
    nodes_faces_vols_2D_coords[:,1:] = vs_reord
    Xs = np.sum(Ns[:,np.newaxis,:,np.newaxis] * nodes_faces_vols_2D_coords[np.newaxis,...],axis=2).transpose(1,0,2)
    #import matplotlib.pyplot as plt
    #plt.scatter(Xs[...,0].flatten(), Xs[...,1].flatten())
    #plt.show
    #plt.savefig("mesh_test.png")

    z_center_abs = np.max(abs(ff[:,:,-1]))/2
    faces_contour_vols_2D_idx_flt = faces_vols[abs(ff[:,:,-1])==z_center_abs]
    xs2 = np.split(faces_contour_vols_2D_idx_flt[...,np.newaxis],ctes.n_volumes)
    faces_around_vols_2D_idx = np.concatenate(xs2,axis=-1).T
    return Xs, nodes_faces_vols_2D_coords, faces_around_vols_2D_idx

def map_SPs_on_faces(Xs, nodes_faces_vols_2D_coords):
    # I think there is an easier way to do that and just for 1 pattern CV:
    # I can do only one for in range 4 or an array filling each row with the matrix pattern I am noticing
    faces_vols_coef = np.empty((ctes.n_volumes,4))
    for i in range(4):
        if i==3: i=-1
        faces_vols_coef[:,i] = (nodes_faces_vols_2D_coords[:,i+1,1]-nodes_faces_vols_2D_coords[:,i,1])/ (nodes_faces_vols_2D_coords[:,i+1,0]-nodes_faces_vols_2D_coords[:,i,0])

    faces_vols_straight = np.zeros_like(faces_vols_coef,dtype=bool)
    faces_vols_straight[(faces_vols_coef==0) + np.isinf(faces_vols_coef)] = True
    SPs_on_faces = np.empty((ctes.n_volumes,n_points,4))

    for i in range(n_points):
        for j in range(4):
            C1 = (Xs[:,i,1] - nodes_faces_vols_2D_coords[:,j,1])
            #C1[abs(C1)<1e-15] = 0
            C2 = faces_vols_coef[:,j]*(Xs[:,i,0] - nodes_faces_vols_2D_coords[:,j,0])
            #C2[abs(C2)<1e-15] = 0
            err = 1e-12
            SPs_on_faces[:,i,j] = abs(C1 - C2)<err
            face_v =  np.isinf(abs(faces_vols_coef[:,j]))
            face_h = (faces_vols_coef[:,j]==0)
            if j==3: j=-1
            #import pdb; pdb.set_trace()

            C_v = abs(Xs[face_v,i,0] - nodes_faces_vols_2D_coords[face_v,j,0])<err
            C_h = abs(Xs[face_h,i,1] - nodes_faces_vols_2D_coords[face_h,j,1])<err
            SPs_on_faces[face_v,i,j] = C_v * \
                (((nodes_faces_vols_2D_coords[face_v,j+1,1]<=Xs[face_v,i,1]) * \
                (Xs[face_v,i,1]<=nodes_faces_vols_2D_coords[face_v,j,1])) + \
                ((nodes_faces_vols_2D_coords[face_v,j+1,1]>=Xs[face_v,i,1]) * \
                (Xs[face_v,i,1]>=nodes_faces_vols_2D_coords[face_v,j,1])) +\
                (abs(Xs[face_v,i,1]-nodes_faces_vols_2D_coords[face_v,j,1])<err) +\
                (abs(Xs[face_v,i,1]-nodes_faces_vols_2D_coords[face_v,j+1,1])<err))
            SPs_on_faces[face_h,i,j] = C_h * \
                (((nodes_faces_vols_2D_coords[face_h,j+1,0]<=Xs[face_h,i,0]) *\
                (Xs[face_h,i,0]<=nodes_faces_vols_2D_coords[face_h,j,0])) + \
                ((nodes_faces_vols_2D_coords[face_h,j+1,0]>=Xs[face_h,i,0]) * \
                (Xs[face_h,i,0]>=nodes_faces_vols_2D_coords[face_h,j,0])) + \
                (abs(Xs[face_h,i,0]-nodes_faces_vols_2D_coords[face_h,j,0])<err) + \
                (abs(Xs[face_h,i,0]-nodes_faces_vols_2D_coords[face_h,j+1,0])<err))
    SPs_on_faces = SPs_on_faces.astype(bool)
    return SPs_on_faces #bool of all SPs

def map_FPs_from_SPs(n_points,faces_contour_vols_2D_idx,SPs_on_faces):
    internal_faces_contour_vols = np.isin(faces_contour_vols_2D_idx,ctes.internal_faces)
    internal_faces_contour_vols_fl = faces_contour_vols_2D_idx[internal_faces_contour_vols]

    st = []
    i=0
    j=0
    idx_2 = np.empty((ctes.n_internal_faces,2),dtype=int)
    for el in internal_faces_contour_vols_fl:
        if el not in st:
            st.append(el)
            idx_2[j,0] = i
            j+=1
        else:
            iel = st.index(el)
            idx_2[iel,1] = i
        i+=1

    v0_SPs_fl = SPs_on_faces.transpose(0,2,1)[internal_faces_contour_vols]

    v0_SPs_2 = v0_SPs_fl[idx_2]
    reorder_internal_faces_to_origin = (np.array(st)[:, None] == ctes.internal_faces).argmax(axis=0)
    v0_SPs = v0_SPs_2[reorder_internal_faces_to_origin]
    return v0_SPs

"""
def get_vols_external_faces(M):
    "function to get the number of external faces of a volume"
    faces_vols = M.volumes.bridge_adjacencies(M.volumes.all,3,2)
    vols_faces_internal_bool = np.isin(faces_vols,ctes.internal_faces)
    vols_n_external_faces = (~vols_faces_internal_bool).sum(axis=-1)
    return vols_n_external_faces
"""
def RT0_shape_functions_2D(n_points,x_points,y_points):
    #function to map from the reference domain to the physical one
    phi = np.empty((n_points,2,2))
    phi[:,0,0] = -1/4*(1-x_points)
    phi[:,0,1] = 1/4*(1+x_points)
    phi[:,1,0] = -1/4*(1-y_points)
    phi[:,1,1] = 1/4*(1+y_points)
    return phi

def GLL(n_points):
    "see if it is going to be necessary. if yes, complete creating the other functions"
    if n_points<2:
        raise ValueError("N points must be larger than 1")
    else:
        x1=-1;
        xn=1;
        w1 = 2/(n*(n-1))
        wn=w1
        x = np.zeros(n_points)
        w = np.zeros(n_points)
        for i in range(n_points):
            x[i] = (1-(3*(n_points-2))/(8*(n_points-1)**3))*cos((4*i-3)/(4*(n_points-1)+1)*math.pi)
            erro_r = 1
            while erro_r<1e-15:
                y = dLgP(n_points-1,x[i])
                y1 = d2LgP(n_points-1,x[i])
                y2 = d3LgP(n_points-1,x[i])
                dx = 2*y*y1/(2*y1**2-y*y2)
                x[i]-=dx
                erro_r=abs(dx)
            w[i] = 2/(n_points*(n_points-1)*lgP(n_points-1,x[i])**2)
        return x, w

def run(M):
    global n_points
    global points
    global weights
    global L
    global dL
    global dgRB
    global dgLB
    global V
    global vols_vec
    global Dr
    global phi
    global mesh_dim
    global Xs

    mesh_dim = 2-1 * (M.volumes.internal.shape==0)

    n_points = data_loaded['compositional_data']['FR']['order']

    #GL = quadpy.c1.gauss_lobatto(n_points)
    if n_points==2:
        points = np.array([-1,1])
        weights = np.array([1,1])

    if n_points==3:
        points = np.array([-1,0,1])
        weights = np.array([1/3, 4/3, 1/3])

    if n_points==4:
        points = np.array([-1,-0.4472136,0.4472136, 1])
        weights = np.array([1/6, 5/6, 5/6, 1/6])

    #points = GL.points
    #weights = GL.weights
    L, dL = Lagrange_poly(n_points, points)
    dgRB, dgLB, gRB, gLB = correction_function(n_points)

    points[abs(points)<1e-15] = 0
    L[abs(L)<5e-15] = 0
    dL[abs(dL)<5e-15] = 0

    #vols_n_external_faces = get_vols_external_faces(M)
    if mesh_dim==2:
        global x_points
        global y_points
        global v0_SPs
        global nFPs
        global reshape_to_points_matrix
        global points_matrix
        global SPs_on_faces

        n_points *= n_points
        n_pmatrx, neta = np.meshgrid(points, points)
        w_pmatrx, weta = np.meshgrid(weights, weights)
        x_points = np.reshape(n_pmatrx,[n_points])
        w_lined = np.reshape(w_pmatrx,[n_points])
        y_points = np.reshape(neta,[n_points])
        w_eta = np.reshape(weta,[n_points])
        points = x_points.copy() #change latter
        Ns = bilinear_shape_functions_2D(n_points,x_points, y_points)
        Xs, nodes_faces_vols_2D_coords, faces_contour_vols_2D_idx = map_faces_and_points_in_vols(M, Ns)
        SPs_on_faces = map_SPs_on_faces(Xs, nodes_faces_vols_2D_coords)
        v0_SPs = map_FPs_from_SPs(n_points,faces_contour_vols_2D_idx,SPs_on_faces)
        phi = RT0_shape_functions_2D(n_points,x_points,y_points)
        nFPs = np.sqrt(n_points).astype(int)
        vols_vec = auxiliary_terms_2D(faces_contour_vols_2D_idx)
        idx = np.linspace(0,n_points-1,n_points,dtype=int)
        aux_reshape = np.split(idx[:,np.newaxis],int(np.sqrt(n_points)),axis=0)
        reshape_to_points_matrix = np.concatenate(aux_reshape,axis=1).T
        points_matrix = x_points[reshape_to_points_matrix]

    else:
        phi = RT0_shape_functions(n_points, points)
        _points, vols_vec = auxiliary_terms(faces_contour_vols_2D_idx)

    V = np.polynomial.legendre.legvander(points,n_points-1)

    #phi = RT0_shape_functions(n_points, points)

    #V_H = Vandermonde2(n_points, points)
    #GV_H = gradVandermonde(n_points, points)
    #Dr = GV_H@np.linalg.inv(V_H)
