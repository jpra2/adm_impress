from packs.manager.mesh_data import MeshData
import os
from packs import defpaths
import scipy.io  as sio
import numpy as np
import matplotlib.pyplot as plt
from packs.manager.meshmanager import create_initial_3D_mesh_prperties, MeshProperty, load_mesh_properties
from packs.examples.load_3d_mesh import first_func, second_func
import pandas as pd
from packs import defpaths
from shutil import rmtree

def load_from_arr(cases: list):
    files = np.array(os.listdir(defpaths.flying))
    df_files = pd.Series(files)
    my_files = []
    for case in cases:
        file = files[(df_files.str.contains(case) & df_files.str.contains('arr'))]        
        my_files.append(os.path.join(defpaths.flying ,file[0]))
    dfs = []
    for f in my_files:
        dfs.append(pd.read_csv(f))
    
    dfs = pd.concat(dfs, ignore_index=True)
    return dfs

def insert_arr_data(df: pd.DataFrame, mesh_data: MeshData, case: str):
    
    remove_loop = ['day', 'case']
    df_case = df[df.case == case]
    day_column = 'day'
    all_days = pd.unique(df_case.day)
    for data_name in df.columns:
        if data_name in remove_loop:
            continue
        data = df_case[data_name].values
        for day in all_days:
            data_day = data[df_case[day_column] == day]
            tag_name = data_name + '_DAY=' + str(day)
            mesh_data.create_tag(tag_name)
            mesh_data.insert_tag_data(tag_name, data_day, 'volumes')

def export_database_from_case(case, dfs, mesh_path):
    df_case = dfs[dfs['case'] == case]
    folder = os.path.join(defpaths.results, case)
    
    if os.path.exists(folder):
        rmtree(folder)
    os.makedirs(folder)
    os.makedirs(os.path.join(folder, 'gifs'))
     
    remove_loop = ['day', 'case']
    
    all_days = pd.unique(df_case['day'])
    for day in all_days:
        
        mesh_data = MeshData(dim=3, mesh_path=mesh_path)        
        df_day = df_case[df_case['day'] == day]
               
        for data_name in df_day.columns:
            if data_name in remove_loop:
                continue
            data = df_day[data_name].values
            mesh_data.create_tag(data_name)
            mesh_data.insert_tag_data(data_name, data, 'volumes')
        
        to_export_name = 'DAY=' + str(day)
        to_export_name = os.path.join(case, to_export_name)
        mesh_data.export_all_elements_type_to_vtk(to_export_name, 'volumes')

def export_data_arr_from_cases1(cases, mesh_properties, mesh_path):
    dfs = load_from_arr(cases)   
    
    for case in cases:
        
        mesh_data = MeshData(dim=3, mesh_path=mesh_path)
        nblocks = len(mesh_properties.volumes)
        ids = np.arange(nblocks)

        mesh_data.create_tag('id', data_type='int')
        mesh_data.insert_tag_data('id', ids, 'volumes')

        permeability = np.ones(nblocks)*100
        mesh_data.create_tag('permeability')
        mesh_data.insert_tag_data('permeability', permeability, 'volumes')

        porosity = np.ones(nblocks)*0.28
        porosity[mesh_properties.volumes_centroids[:, 0] > 720] = 0.3
        mesh_data.create_tag('porosity')
        mesh_data.insert_tag_data('porosity', porosity, 'volumes')
        
        to_export_name = case
        insert_arr_data(dfs, mesh_data, case)
        mesh_data.export_all_elements_type_to_vtk(to_export_name, 'volumes')
        export_database_from_case(case, dfs, mesh_path)


def create_gif(case):
    
    import glob
    import contextlib
    from PIL import Image
    
    gif_path = os.path.join(defpaths.results, case)
    gif_path = os.path.join(gif_path, 'gifs')

    # filepaths
    fp_in = gif_path
    fp_out = os.path.join(gif_path, 'movie.gif')

    # use exit stack to automatically close opened images
    with contextlib.ExitStack() as stack:
        all_dirs = os.listdir(fp_in)
        all_dirs = [os.path.join(fp_in, i) for i in all_dirs]

        # lazily load images
        imgs = (stack.enter_context(Image.open(f))
                for f in sorted(all_dirs))

        # extract  first image from iterator
        img = next(imgs)

        # https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
        img.save(fp=fp_out, format='GIF', append_images=imgs,
                save_all=True, duration=200, loop=0)
    
    
    # files = os.listdir(gif_path)
    # images = []
    # for filename in files:
    #     file_img = os.path.join(gif_path, filename)
    #     images.append(imageio.imread(file_img))
    
    # gif_path = os.path.join(gif_path, 'movie.gif')
    # imageio.mimsave(gif_path, images)
            

    

# #######################################################
# mesh_name = '9x9x1_sim.h5m'
# mesh_path = os.path.join(defpaths.mesh, mesh_name)
mesh_properties_name = '9x9x1_sim'

# first_func(mesh_path, mesh_properties_name)
# second_func(mesh_properties_name)
mesh_properties = load_mesh_properties(mesh_properties_name)

# cases = ['base12', 'base13', 'base14', 'base15', 'base16', 'base17']
# export_data_arr_from_cases1(cases, mesh_properties, mesh_path)

case = 'base12'
create_gif(case)
print(mesh_properties)



