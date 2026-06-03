import os
import pandas as pd

flying = 'flying'
data_mesh = os.path.join(flying, 'data_mesh_property')
data_simulation = data_mesh

plots_folder = os.path.join(flying, 'plots')
results = 'results'
# pressure_results = os.path.join(results, 'pressure_ameba')
# saturation_results = os.path.join(results, 'saturation_ameba')
# pressure_results_ms = os.path.join(results, 'pressure_ameba_ms')
# saturation_results_ms = os.path.join(results, 'saturation_ameba_ms')

pressure_results = os.path.join(results, 'pressure_sin')
saturation_results = os.path.join(results, 'saturation_sin')
pressure_results_ms = os.path.join(results, 'pressure_sin_ms')
saturation_results_ms = os.path.join(results, 'saturation_sin_ms')

# pressure_results = os.path.join(results, 'pressure_het')
# saturation_results = os.path.join(results, 'saturation_het')
# pressure_results_ms = os.path.join(results, 'pressure_het_ms')
# saturation_results_ms = os.path.join(results, 'saturation_het_ms')

# pressure_results = os.path.join(results, 'pressure_layers')
# saturation_results = os.path.join(results, 'saturation_layers')
# pressure_results_ms = os.path.join(results, 'pressure_layers_ms')
# saturation_results_ms = os.path.join(results, 'saturation_layers_ms')
# pressure_results_ms = os.path.join(results, 'pressure_layers_ms_amsu')
# saturation_results_ms = os.path.join(results, 'saturation_layers_ms_amsu')

mesh = 'mesh'
data_folder = 'data'
mpfad_mesh_folder = 'mpfad_mesh_tests'
lpew2_mesh_folder = 'lpew2_mesh_test'
unstructured_coarse_test_mesh_folder = 'uns_coarse_test'

remove_folder = 'remove'

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

matrix_path = os.path.join(flying, 'matrices.h5')

layers_folder = os.path.join(unstructured_coarse_test_mesh_folder, 'layers')
finescale_mesh_layer = os.path.join(layers_folder, '5spotHETq.msh')
coarse_mesh_layer = os.path.join(layers_folder, 'coarse1.msh')
layers_tri_finescale = os.path.join(layers_folder, 'layers_tri_f.msh')
layers_tri_coarse = os.path.join(layers_folder, 'layers_tri_coarse.msh')
layers_tri_finescale_ff = os.path.join(layers_folder, 'layers_tri_ff.msh')
layers_tri_coarse_f = os.path.join(layers_folder, 'layers_tri_coarse_f.msh')
coarse1_layers = os.path.join(layers_folder, 'coarse_test.msh')
fine1_layers = os.path.join(layers_folder, 'fine1_layers.msh')
coarse3_layers = os.path.join(layers_folder, 'coarse3.msh')
fine3_layers = os.path.join(layers_folder, 'fine3_layers.msh')
coarse4_layers = os.path.join(layers_folder, 'coarse4.msh')
fine4_layers = os.path.join(layers_folder, 'fine4.msh')
coarse5_layers = os.path.join(layers_folder, 'coarse5.msh')
fine5_layers = os.path.join(layers_folder, 'fine5.msh')

ameba_folder = os.path.join(unstructured_coarse_test_mesh_folder, 'ameba')
ameba_fine1 = os.path.join(ameba_folder, 'amebaquad1.msh')
ameba_fine2 = os.path.join(ameba_folder, 'amebafine2.msh')
ameba_fine2_v4 = os.path.join(ameba_folder, 'amebafine2_v4.msh')
ameba_coarse2 = os.path.join(ameba_folder, 'coarse2_2_ameba.msh')
ameba_fine3 = os.path.join(ameba_folder, 'amebafine3.msh')
ameba_fine4 = os.path.join(ameba_folder, 'amebafine4.msh')
ameba_coarse4f = os.path.join(ameba_folder, 'coarse4f_ameba.msh')
ameba_fine5 = os.path.join(ameba_folder, 'amebafine5.msh')


biphasic_folder = 'bifasico'
barreira_mesh = os.path.join(biphasic_folder, '5spotN2f.msh')
barreira_mesh_coarse = os.path.join(biphasic_folder, '5spotN2_coarse_2.msh')
quadrado_mesh = os.path.join(biphasic_folder, '5spot_quadrado_no_meio_v2.msh')
quadrado_mesh_coarse = os.path.join(biphasic_folder, '5spot_quadrado_no_meio_v2_coarse.msh')
ex1_finescale = os.path.join(biphasic_folder, 'ex1f.msh')
sin1_fine = os.path.join(biphasic_folder, 'sin1_fine.msh')
sin1_fine_v4 = os.path.join(biphasic_folder, 'sin1_fine_v4.msh')
sin1_coarse = os.path.join(biphasic_folder, 'sin1_coarse.msh')
sin1_coarse2 = os.path.join(biphasic_folder, 'sin1_coarse2.msh')

het_fine1_mesh = os.path.join(biphasic_folder, 'het_fine1.msh')
het_coarse1_mesh = os.path.join(biphasic_folder, 'het_coarse1.msh')
het_fine2_mesh = os.path.join(biphasic_folder, 'het_fine2.msh')
het_coarse2_mesh = os.path.join(biphasic_folder, 'het_coarse2.msh')

symetric_finescale = os.path.join(biphasic_folder, 'symetricf.msh')
symetric_finescale_tri = os.path.join(biphasic_folder, 'symetricf_tri.msh')
symetric_coarse = os.path.join(biphasic_folder, 'symetric_coarse.msh')

points_chueh_artur_path = os.path.join(data_folder, 'points_chueh_artur.mat')

ftest1 = 'malha_teste_1.msh'
ftest1_v4 = 'malha_teste_1_v4.msh'
ctest1_1 = 'malha_teste_1_coarse_1.msh'
ctest1_2 = 'malha_teste_1_coarse_2.msh'
ctest1_3 = 'malha_teste_1_coarse_3.msh'
ctest1_4 = 'malha_teste_1_coarse_4.msh'
ctest1_5 = 'malha_teste_1_coarse_5.msh'

ftest3 = 'malha_teste_3.msh'
ftest3_v4 = 'malha_teste_3_v4.msh'