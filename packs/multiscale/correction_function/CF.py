from scipy.sparse import csc_matrix
import numpy as np

def solve_correction(LU,q):

    gids=[]
    cf=[]
    for lu, gid in zip(LU[0],LU[1]):
        gids.append(gid)
        cf.append(lu.solve(q[gid]))


    gids=np.concatenate(gids)
    cf=np.concatenate(cf)
    cf_map = np.zeros_like(q,dtype=float)
    cf_map[gids] = cf
    return sorted(gids), cf_map

def get_cf_f(LUf, q, cfe, gidse, Tfe):
    gidsf, cf_f_f=solve_correction(LUf,q)
    qe=np.zeros_like(cfe)
    qe[gidsf]=-Tfe*cfe[gidse]
    gidsf, cf_f_e=solve_correction(LUf,qe)
    return gidsf, cf_f_f, cf_f_e

def get_cf_i(LUi, q, cf_e_e, cf_f_f, cf_f_e, gids_e, gids_f, Tif, Tie):
    gidsi, cf_i_i=solve_correction(LUi,q)
    qf=np.zeros_like(cf_f_f)
    qf[gidsi]=-Tif*cf_f_f[gidsf]
    gidsi, cf_i_f=solve_correction(LUi,qf)
    qe = np.zeros_like(cf_f_f)
    qe[gids_i]=-Tie*cf_e_e[gids_e]-Tif*cf_f_e[gids_e]
    gidsi, cf_i_e=solve_correction(LUi,qe)
    return cf_i_i+cf_i_f+cf_i_e

def get_correction_function(LocalLU, cMats, q):
    # assert len(cMats) == 6
    Tfe = cMats['Afe']
    Tif = cMats['Aif']
    Tie = cMats['Aie']

    gids_e, cf_e_e = solve_correction(LocalLU[2], q)
    gids_f, cf_f_f, cf_f_e = get_cf_f(LocalLU[1], q, cf_e_e, gids_e, Tfe)
    if len(LocalLU[0])>0:
        cf_i = get_cf_i(LocalLU[0], q, cf_e_e, cf_f_f, cf_f_e, gids_e, gids_f, Tif, Tie)
    else:
        cf_i = np.zeros_like(cf_e_e)

    return cf_i + cf_f_f + cf_f_e + cf_e_e
