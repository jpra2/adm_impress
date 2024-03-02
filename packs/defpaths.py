import os
import pandas as pd

flying = 'flying'
results = 'results'
mesh = 'mesh'
mpfad_mesh_folder = 'mpfad_mesh_tests'
lpew2_mesh_folder = 'lpew2_mesh_test'
unstructured_coarse_test_mesh_folder = 'uns_coarse_test'

def load_mpfad_meshs_by_name(mesh_name: str, folder: str) -> pd.DataFrame:
    file_names = os.listdir(folder)
    names_df = pd.Series(data=file_names)
    my_meshs = names_df[names_df.str.contains(mesh_name)]
    my_meshs = [os.path.join(mpfad_mesh_folder, name) for name in my_meshs.values]
    mesh_n = [mesh_name]*len(my_meshs)
    df_mesh = pd.DataFrame({
        'mesh_path':  my_meshs,
        'mesh_type': mesh_n
    })
    
    return df_mesh
        
def load_mpfad_meshs() -> pd.DataFrame:
    global mesh, mpfad_mesh_folder
    folder = os.path.join(mesh, mpfad_mesh_folder)
    meshs_header = ['mesh1', 'mesh2', 'mesh5', 'mesh6']
    dfs = []
    for mesh_name in meshs_header:
        dfs.append(
            load_mpfad_meshs_by_name(mesh_name, folder)
        )
    
    dfs = pd.concat(dfs, ignore_index=True)
    return dfs
    
def load_mpfad_meshtest_by_type_and_number(mesh_type: str, n: int):
    meshs_df = load_mpfad_meshs()
    meshs_df = meshs_df[meshs_df['mesh_type'] == mesh_type]
    meshs_df = meshs_df[meshs_df['mesh_path'].str.contains('_' + str(n))]
    mesh_path = meshs_df['mesh_path'].values[0]
    return mesh_path
    
def load_su_mesh_paths():
    mesh_test_su_mpfa = 'mesh_su'
    global mesh, mpfad_mesh_folder
    folder = os.path.join(mesh, mpfad_mesh_folder)
    file_names = os.listdir(folder)
    names_df = pd.Series(data=file_names)
    su_meshs = names_df[names_df.str.contains(mesh_test_su_mpfa)].values
    su_meshs = [os.path.join(mpfad_mesh_folder, mesh_path) for mesh_path in su_meshs]
    return su_meshs
    
    
    
mpfad_test_mesh = '2d_unstructured.msh'
mpfad_mesh_properties_name = 'gls_test_weights'
mpfad_mesh_properties_neumann_name = 'neumann_gls_test_weights'

mpfad_mesh_2d_test_6 = '2d_test6_paper.h5m'
mesh_properties_2d_test_6_name = 'test6_mpfad'

oblique_quad_mesh = os.path.join(mpfad_mesh_folder, 'oblique_quadrilateral_test1.msh')
mesh_prop_test7 = 'test_oblique_7_mpfad'
mesh_properties_mesh1_test = 'mpfad_lsds_mesh1'

linear_2k_test = os.path.join(mpfad_mesh_folder, 'linear_2k_test.msh')
mesh_prop_linear_2k = 'linear_2k_test'


