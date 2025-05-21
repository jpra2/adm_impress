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
     















