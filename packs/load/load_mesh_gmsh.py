import gmsh
import numpy as np
import os

def run():

    file_name = os.path.join('mesh', 'uns_coarse_test', 'mesh0_3.msh')
    gmsh.initialize()
    gmsh.open(file_name)

    dim=-1
    tag = -1
    nodeTags, coords, parametricCoords = gmsh.model.mesh.getNodes(dim=dim, tag=tag)

    

    # gmsh.fltk.run()
    gmsh.finalize()