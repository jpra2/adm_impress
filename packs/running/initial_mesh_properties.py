
from ..data_class.data_impress import Data
from ..data_class.elements_lv0 import ElementsLv0
from ..contours.wells import Wells
from ..convert_unit.conversion import Conversion
from ..preprocess.preprocess1 import set_saturation_regions
from ..preprocess.prep0_0 import Preprocess0
from ..directories import data_loaded, only_mesh_name, name_variable_inputs_file_load
from ..directories import simulation_type, types_simulation
from ..multiscale.preprocess.dual_primal.create_dual_and_primal_mesh import MultilevelData
import numpy as np
import os
from .update_inputs import Update_variable, Delete_Files
# import time
import pdb

# import pdb; pdb.set_trace()

def initial_mesh():

    global data_loaded

    load = data_loaded['load_data']
    convert = data_loaded['convert_english_to_SI']
    multilevel_data = data_loaded['multilevel_data']
    load_multilevel_data = data_loaded['load_multilevel_data']
    Update_variable(name_variable_inputs_file_load)
    Delete_Files()

    if multilevel_data and load_multilevel_data:
        import time
        t0=time.time()
        print("creating M")
        from ..load.preprocessor_load import init_mesh
        # M = init_mesh('flying/multilevel_data-all.h5m')
        import time

        M = init_mesh('saves/initial_mesh_' + only_mesh_name + '.h5m')
        print("time to create M: {} seconds".format(time.time()-t0))
        t0=time.time()
        elements_lv0 = ElementsLv0(M, load=load)
        print("time_to_create elements_lv0: {}".format(time.time()-t0))
        print("creating data impress")
        data_impress = Data(M, elements_lv0, load=load)
        print("time_to_create data_impress: {}".format(time.time()-t0))
        print("creating ml_data")
        t0=time.time()
        if not load:
            print("preprocess0")
            t0=time.time()
            Preprocess0(M, elements_lv0)
            print(time.time()-t0,"preprocess0")
        ml_data = MultilevelData(data_impress, M, load=load_multilevel_data)
        ml_data.load_tags()
        print("time_to_create ml_data: {}".format(time.time()-t0))

    else:
        if load:
            import time
            print("creating M")
            t0=time.time()
            from ..load.preprocessor_load import init_mesh
            M = init_mesh('saves/initial_mesh_' + only_mesh_name + '.h5m')
            print("time to create M: {} seconds".format(time.time()-t0))
        else:
            import time
            print("creating M")
            t0=time.time()
            from ..load.preprocessor0 import M
            print("time to create M: {} seconds".format(time.time()-t0))
        t0=time.time()
        elements_lv0 = ElementsLv0(M, load=load)
        print("time_to_create elements_lv0: {}".format(time.time()-t0))
        t0=time.time()
        data_impress = Data(M, elements_lv0, load=load)
        print("time_to_create elements_data_impress: {}".format(time.time()-t0))
        t0=time.time()
        if not load:
            Preprocess0(M, elements_lv0)
        if multilevel_data:
            ml_data = MultilevelData(data_impress, M)
            ml_data.run()
            # ml_data.save_mesh()
        print("time_to_create ml_data: {}".format(time.time()-t0))

    if simulation_type not in types_simulation:
        raise Error('Invalid simulation type')

    if simulation_type == 'compositional':
        from ..contours.wells import WellsCompositional
        wells = WellsCompositional(M, load = load)
    else:
        wells = Wells(M, elements_lv0, load=load)

    biphasic = data_loaded['biphasic']
    load_biphasic_data = data_loaded['load_biphasic_data']

    if biphasic and not load_biphasic_data:
        set_saturation_regions(M, wells)
        # TODO: atualizar gama

    if convert and not load:
        conversion = Conversion(wells, data_impress)
        conversion.convert_English_to_SI()
        mi_monophasic = data_loaded['monophasic_data']['mi']
        data_impress['transmissibility'] = data_impress['pretransmissibility']*np.repeat(mi_monophasic, len(data_impress['pretransmissibility']))

    if not load:

        wells.update_values_to_mesh()
        data_impress.update_variables_to_mesh()
        wells.export_all_datas_to_npz()
        M.save_variables('initial_mesh_'+ only_mesh_name)
        del data_impress['permeability']

    if data_loaded['deletar_results'] == True:
        results_file = 'results'
        for arq in os.listdir(results_file):
            if arq.endswith('.vtk'):
                os.remove(os.path.join(results_file, arq))

    return M, elements_lv0, data_impress, wells
