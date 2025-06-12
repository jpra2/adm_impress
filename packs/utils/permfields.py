import numpy as np
from numpy.random import RandomState
from packs.utils.utils_old import is_point_inside_circle
from packs import defpaths
import scipy.io as sio

def _chueh_2(centroid, N_centroids, dists, phiL):
        
        dists[:] = np.linalg.norm(centroid - N_centroids, axis=1)
        phiL[:] = np.exp(-1*np.power(dists/0.05, 2))
        v1 = max([phiL.sum(), 0.01])
        v2 = min([v1, 4.0])
        return v2

def random_permeability_chueh(elements_centroids: np.ndarray, N: int, state: int = 2, aditional_ids=np.array([])):
    gen = RandomState(state)
    nelements = elements_centroids.shape[0]
    random_choice = gen.randint(0, nelements, size=N)
    if aditional_ids.shape[0] > 0:
        random_choice = np.union1d(random_choice, aditional_ids)

    r1 = 0.07;
    r2 = 0.07;
    x1 = 0.3;
    y1 = 0.3;
    x2 = 1.7;
    y2 = 0.5;

    ids_to_remove = []
    for pid in random_choice:
         pcentroid = elements_centroids[pid]
         is_in_circle1 = is_point_inside_circle(
              pcentroid[0],
              pcentroid[1],
              x1,
              y1,
              r1
         )
         is_in_circle2 = is_point_inside_circle(
              pcentroid[0],
              pcentroid[1],
              x2,
              y2,
              r2
         )

         if is_in_circle1 == True or is_in_circle2 == True:
              ids_to_remove.append(pid)
    
    ids_to_remove = np.array(ids_to_remove)


    
    if ids_to_remove.shape[0] > 0:
         random_choice = np.setdiff1d(random_choice, ids_to_remove)

    perm = np.zeros(nelements)
    N_centroids = elements_centroids[random_choice]
    dists = np.zeros(random_choice.shape[0])
    phiL = dists.copy()

    for i in range(nelements):
        perm[i] = _chueh_2(elements_centroids[i], N_centroids, dists, phiL)
    
    return perm

def random_permeability_chueh_v2(elements_centroids: np.ndarray, N: int, state: int = 2, aditional_ids=np.array([])):
    gen = RandomState(state)
    nelements = elements_centroids.shape[0]
    random_choice = gen.randint(0, nelements, size=N)
    if aditional_ids.shape[0] > 0:
        random_choice = np.union1d(random_choice, aditional_ids)

    r1 = 0.07
    r2 = 0.07
    r3 = 0.07
    r4 = 0.07
    r5 = 0.07
    r6 = 0.07
    r7 = 0.07
    x1 = 0.3
    y1 = 0.3
    x2 = 1.7
    y2 = 0.5
    x3 = 1.38
    y3 = 0.395
    x4 = 1.36
    y4 = 0.825
    x5 = 1.245
    y5 = 0.855
    x6 = 1.735
    y6 = 0.63
    x7 = 1.624
    y7 = 0.454

    remove = False

    x_list = [x1, x2, x3, x4]
    y_list = [y1, y2, y3, y4]
    r_list = [r1, r2, r3, r4]



    ids_to_remove = []

    for pid in random_choice:
         
         pcentroid = elements_centroids[pid]

         is_in_circle1 = is_point_inside_circle(
              pcentroid[0],
              pcentroid[1],
              x1,
              y1,
              r1
         )
         is_in_circle2 = is_point_inside_circle(
              pcentroid[0],
              pcentroid[1],
              x2,
              y2,
              r2
         )
         is_in_circle3 = is_point_inside_circle(
                    pcentroid[0],
                    pcentroid[1],
                    x3,
                    y3,
                    r3
         )

         is_in_circle4 = is_point_inside_circle(
                    pcentroid[0],
                    pcentroid[1],
                    x4,
                    y4,
                    r4
         )

         is_in_circle5 = is_point_inside_circle(
                    pcentroid[0],
                    pcentroid[1],
                    x5,
                    y5,
                    r5
         )

         is_in_circle6 = is_point_inside_circle(
                    pcentroid[0],
                    pcentroid[1],
                    x6,
                    y6,
                    r6
         )

         is_in_circle7 = is_point_inside_circle(
                    pcentroid[0],
                    pcentroid[1],
                    x7,
                    y7,
                    r7
         )


         verify = is_in_circle1 or is_in_circle2 or is_in_circle3 or is_in_circle4 or is_in_circle5 or is_in_circle6 or is_in_circle7

         if verify == True:
              ids_to_remove.append(pid)
    
    ids_to_remove = np.array(ids_to_remove)


    
    if ids_to_remove.shape[0] > 0:
         random_choice = np.setdiff1d(random_choice, ids_to_remove)

    perm = np.zeros(nelements)
    N_centroids = elements_centroids[random_choice]
    dists = np.zeros(random_choice.shape[0])
    phiL = dists.copy()

    for i in range(nelements):
        perm[i] = _chueh_2(elements_centroids[i], N_centroids, dists, phiL)
    
    return perm

def chueh_perm_artur_paper(elements_centroids: np.ndarray):
     points_def = sio.loadmat(defpaths.points_chueh_artur_path)['points'][:, 1:3]
     nelements = elements_centroids.shape[0]
     perm = np.zeros(nelements)
     
     N_centroids = points_def
     dists = np.zeros(N_centroids.shape[0])
     phiL = dists.copy()

     for i in range(nelements):
        perm[i] = _chueh_2(elements_centroids[i], N_centroids, dists, phiL)
    
     return perm
     
def chueh_1(elements_centroids: np.ndarray):
     xpoints = elements_centroids[:, 0]
     ypoints = elements_centroids[:, 1]

     v1 = (ypoints -0.5 - 0.1*np.sin(10*xpoints))/0.1
     v2 = -1*np.power(v1, 2)
     v3 = np.exp(v2)
     v3[v3 < 0.01] = 0.01

     return v3















