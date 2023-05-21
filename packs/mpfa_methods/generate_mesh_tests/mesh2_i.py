import numpy as np
from pymoab import core, types, rng, topo_util
from packs import defpaths
from packs.utils.utils_old import get_box
import os


def generate_first_points(dxy):
    points = np.array([
        [0,0],
        [dxy, 0],
        [dxy, dxy],
        [0, dxy]
    ])
    
    n_diag1 = 3
    n_diag2 = 4
    
    v_diag2 = points[2] - points[0] # principal diagonal
    v_diag1 = points[1] - points[3] # sec diagonal
    
    vv_diag1 = v_diag1/n_diag1
    vv_diag2 = v_diag2/n_diag2
    
    points_along_y0 = np.array([
        [0, 0],
        [dxy/2, 0],
        [dxy, 0]
    ])
    
    points_along_y1 = np.array([
        [0, dxy],
        [dxy/2, dxy],
        [dxy, dxy]
    ])
    
    points_along_x0 = np.array([
        [0,0],
        [0, dxy/2],
        [0, dxy]
    ])
    
    points_along_x1 = np.array([
        [dxy, 0],
        [dxy, dxy/2],
        [dxy, dxy]
    ])
    
    points_diag1 = np.array([
        [0, dxy],
        [0, dxy] + vv_diag1,
        [0, dxy] + 2*vv_diag1,
        [0, dxy] + 3*vv_diag1
    ])
    
    points_diag2 = np.array([
        [0, 0],
        [0, 0] + vv_diag2,
        [0, 0] + 3*vv_diag2,
        [0, 0] + 4*vv_diag2
    ])
    
    points = np.array([
        points_along_x0[0], #0
        points_along_y0[1], #1
        points_along_y0[2], #2
        points_along_x1[1], #3
        points_along_x1[2], #4
        points_along_y1[1], #5
        points_along_y1[0], #6
        points_along_x0[1], #7
        points_diag2[1],    #8
        points_diag2[2],    #9
        points_diag1[1],    #10
        points_diag1[2]     #11
    ])
    
    return points

def create_tris0(points):
    
    z_points = np.zeros(len(points)).reshape((1, len(points)))
    points2 = np.hstack([points, z_points.T])
    index_points = [
        [0, 1, 8],
        [8, 1, 11],
        [11, 1, 2],
        [11, 2, 3],
        [11, 3, 9],
        [9, 3, 4],
        [9, 4, 5],
        [9, 5, 10],
        [10, 5, 6],
        [10, 6, 7],
        [10, 7, 8],
        [8, 0, 7],
        [8, 11, 10],
        [10, 11, 9]    
    ]
    
    tris0 = points2[index_points]
    
    return tris0

def create_squares(tris0, dxy, n):
    
    delta = dxy/12
    mesh_points0 = tris0.reshape((tris0.shape[0]*tris0.shape[1], -1))
    
    continue_for = set()
    for i, point in enumerate(mesh_points0):
        if i in continue_for:
            continue
        limites = np.array([
            point - delta,
            point + delta
        ])
        index = get_box(mesh_points0, limites)
        if len(index) > 1:
            [continue_for.add(j) for j in index[1:]]
    
    to_remove = np.full(mesh_points0.shape[0], False, dtype=bool)
    index_to_remove = np.array(list(continue_for))
    to_remove[index_to_remove] = True
    to_maintain = ~to_remove
    mesh_points0 = mesh_points0[to_maintain]
    
    vy = np.array([0, dxy, 0])
    vx = np.array([dxy, 0, 0])
    
    points_to_replicate_x = mesh_points0[mesh_points0[:, 0] != 0]
    
    all_mesh_points = [mesh_points0]
    
    [all_mesh_points.append(points_to_replicate_x + i*vx) for i in range(1, n)]
    all_mesh_points = np.concatenate(all_mesh_points)
    
    points_to_replicate_y = all_mesh_points[all_mesh_points[:, 1]!=0]
    
    other_points = []
    for i in range(1, n):
        other_points.append(points_to_replicate_y + i*vy)
    
    other_points = np.concatenate(other_points)
    all_mesh_points = np.vstack([all_mesh_points, other_points])    
    mesh_points_index = np.arange(len(all_mesh_points))
    
    all_tris = []
    
    for j in range(n):
        for i in range(n):
            all_tris.append(tris0 + i*vx + j*vy)
    
    all_tris = np.array(all_tris)
    all_tris = all_tris.reshape((all_tris.shape[0]*all_tris.shape[1], -1, 3))
    
    all_tri_points = []
    for tri in all_tris:
        points_mesh = []
        for point in tri:
            dists = np.linalg.norm(point - all_mesh_points, axis=1)
            point_mesh = mesh_points_index[dists <= dists.min()][0]
            points_mesh.append(point_mesh)
        points_mesh = np.array(points_mesh)
        all_tri_points.append(points_mesh)
    
    all_tri_points = np.array(all_tri_points)
            
    
    return mesh_points_index, all_mesh_points, all_tri_points 
    
def create_mesh(mesh_points_index, all_mesh_points, all_tri_points):
    
    mb = core.Core()
    verts = mb.create_vertices(all_mesh_points.flatten())
    all_tris = []
    all_tri_points = all_tri_points.astype(int)
    for tripoints in all_tri_points:
        verts_tri = [verts[int(i)] for i in tripoints]
        all_tris.append(verts_tri)
    mb.create_elements(types.MBTRI, all_tris) 
    mesh_name = os.path.join(defpaths.mesh, defpaths.mpfad_mesh_2d_test_6)   
    mb.write_file(mesh_name)
    
 

def generate(n=2):
    """This mesh is found in paper:
        Benchmark on Discretization Schemes for Anisotropic Diffusion Problems on General Grids 
        n: number of squares in x and y directions
    """
    
    L = 1
    dxy = L/n
    tris0 = create_tris0(
        generate_first_points(dxy)
    )
    mesh_points_index, all_mesh_points, all_tri_points = create_squares(tris0, dxy, n)
    create_mesh(mesh_points_index, all_mesh_points, all_tri_points)
    
    
    
    
    