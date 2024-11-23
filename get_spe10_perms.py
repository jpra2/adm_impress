import numpy as np
from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh

M = msh('mesh/60x220x1.h5m')

data = np.load('spe10_perms_and_phi.npz')
perms = data['perms']
kx = perms[:,0]
c85 = kx[60*220*84:] # camada 85
c1 = kx[0:60*220] # camada 1

centroids = np.zeros((60*220, 3))
cents = M.volumes.center[:]


for j in range(220):
    for i in range(60):
        centroids[i + 60*j] = [20/2 + 20*i, 10/2 + 10*j, 1]

for i, ci in enumerate(centroids):
    test = (cents[:,0] == ci[0]) & (cents[:,1] == ci[1]) & (cents[:,2] == ci[2])
    M.poro[test] = c85[i]

meshset = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset, M.core.all_volumes)
M.core.mb.write_file('testtt.vtk', [meshset])
    