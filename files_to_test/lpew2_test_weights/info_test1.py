## test weights test from matlab code Fernando
from packs import defpaths
import os

def get_filenames():

    my_path_file = os.path.dirname(os.path.relpath(__file__))
    
    mesh_name = os.path.join(defpaths.lpew2_mesh_folder, 'BenchHydraulic2.msh')
    mesh_properties_name = 'bench_hidraulic'
    matfile = os.path.join(my_path_file, 'resp.mat')

    return {
        'mesh_name': mesh_name,
        'mesh_properties_name':  mesh_properties_name,
        'matfile': matfile
    }