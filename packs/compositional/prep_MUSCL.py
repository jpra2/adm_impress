from packs.utils import constants as ctes
import scipy.sparse as sp
import numpy as np



def get_neiboring_properties_to_reconstruction(M):
    neig_vols = M.volumes.bridge_adjacencies(M.volumes.all,2,3)
    
    lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
    cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
    data = np.array([np.ones(len(ctes.v0[:, 0])), np.ones(len(ctes.v0[:, 0])),
                    np.zeros(len(ctes.v0[:, 0])), np.zeros(len(ctes.v0[:, 0]))]).flatten()
    all_neig = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)).toarray()
    all_neig = all_neig.astype(int)
    all_neig2 = all_neig + np.identity(ctes.n_volumes)
    allneig_and_vol = all_neig2.astype(int)

    pos_neig = M.data['centroid_volumes'].T[:,np.newaxis,:] * allneig_and_vol[np.newaxis,:,:]

    pos = pos_neig.transpose(0,2,1)

    ds = pos_neig - pos
    ds_norm = np.linalg.norm(ds, axis=0)
    versor_ds = np.empty(ds.shape)
    versor_ds[:,ds_norm==0] = 0
    versor_ds[:,ds_norm!=0] = ds[:,ds_norm!=0] / ds_norm[ds_norm!=0]

    ds_vols = ds * versor_ds
    ds_vols = ds_vols.sum(axis = 2)

    all_neig_by_axes = np.repeat(all_neig[np.newaxis,:,:],3,axis=0)
    all_neig_by_axes*=abs(versor_ds).astype(int)
    all_neig_by_axes = all_neig_by_axes.sum(axis=1)

    all_neig = all_neig.sum(axis=1)

    return all_neig, allneig_and_vol, all_neig_by_axes, versor_ds, ds_vols

def identify_contour_faces(all_neig):
    contour_neig = np.min(all_neig)
    vols_contour = np.argwhere(all_neig==contour_neig).flatten()
    faces_contour = []

    vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
    lines = np.arange(ctes.n_internal_faces)
    vols_vec[ctes.v0[:,0],1] = lines
    vols_vec[ctes.v0[:,1],0] = lines

    for i in range(len(vols_contour)):
        f0 = np.argwhere(ctes.v0[:,0] == vols_contour[i]).flatten()
        f1 = np.argwhere(ctes.v0[:,1] == vols_contour[i]).flatten()
        faces_contour.extend(f0.tolist())
        faces_contour.extend(f1.tolist())
    faces_contour = np.array(faces_contour)
    return faces_contour

def run(M):
    global all_neig
    global allneig_and_vol
    global all_neig_by_axes
    global versor_ds
    global ds_vols
    global faces_contour
    global ds_face
    global versor_ds_face
    global ds_face_abs

    ds_face = (M.data['centroid_volumes'][ctes.v0[:,1],:] -  M.data['centroid_volumes'][ctes.v0[:,0],:])
    #if any(ds_face.flatten()<0): import pdb; pdb.set_trace()
    versor_ds_face = np.sign(ds_face.sum(axis=-1))
    ds_face_abs = abs(ds_face)
    all_neig, allneig_and_vol, all_neig_by_axes, versor_ds, ds_vols = \
        get_neiboring_properties_to_reconstruction(M)
    faces_contour = identify_contour_faces(all_neig)
