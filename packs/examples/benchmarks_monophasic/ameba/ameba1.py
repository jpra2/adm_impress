from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths
from packs.examples.benchmarks_monophasic.cross.test_cross_1 import(
    create_primal_ids,
    export_primal_ids,
    create_dual_ids,
    export_dual_ids
)


def get_properties():
    
    fine_mesh_path = defpaths.ameba_fine2
    fine_mesh_properties_name = 'ameba_fine2' 
    fine_mesh_path_v4 = defpaths.ameba_fine2_v4

    coarse_mesh_path = defpaths.ameba_coarse2
    coarse_mesh_properties_name = 'ameba_coarse2'

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name, mesh_name_v4=fine_mesh_path_v4)
    coarse_properties = preprocess_mesh(coarse_mesh_path, coarse_mesh_properties_name)

    return fine_properties, coarse_properties, fine_mesh_path, coarse_mesh_path


def run4():

    fp, cp, fine_mesh_path, coarse_mesh_path = get_properties()
    
    bool_export_primal_id = True
    bool_export_dual_id = True
    my_dual_type = 1

    create_primal_ids(fp, cp, update=bool_export_primal_id)
    export_primal_ids(fine_mesh_path, fp, coarse_mesh_path, export=bool_export_primal_id)
    create_dual_ids(fp, cp, update=bool_export_dual_id, dual_type=my_dual_type)
    export_dual_ids(fine_mesh_path, fp, export=bool_export_dual_id)


    import pdb; pdb.set_trace()
