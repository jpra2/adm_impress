
from ..data_class.data_impress import Data
from ..data_class.elements_lv0 import ElementsLv0
from ..contours.wells import Wells
from ..convert_unit.conversion import Conversion
from ..preprocess.preprocess1 import set_saturation_regions
from ..preprocess.prep0_0 import Preprocess0
from ..directories import data_loaded
from ..multiscale.preprocess.dual_primal.create_dual_and_primal_mesh import MultilevelData
import numpy as np

# import pdb; pdb.set_trace()

def initial_mesh(load=False, convert=False):

    global data_loaded

    multilevel_data = data_loaded['multilevel_data']
    load_multilevel_data = data_loaded['load_multilevel_data']

    if multilevel_data and load_multilevel_data:
        from ..load.preprocessor_load import init_mesh
        M = init_mesh('flying/multilevel_data-all_.h5m')
        elements_lv0 = ElementsLv0(M, load=load)
        data_impress = Data(M, elements_lv0, load=load)
        if not load:
            Preprocess0(M, elements_lv0)
        ml_data = MultilevelData(data_impress, M, load=load_multilevel_data)
        ml_data.load_tags()

    else:
        from ..load.preprocessor0 import M
        elements_lv0 = ElementsLv0(M, load=load)
        data_impress = Data(M, elements_lv0, load=load)
        if not load:
            Preprocess0(M, elements_lv0)
        if multilevel_data:
            ml_data = MultilevelData(data_impress, M)
            ml_data.run()
            ml_data.save_mesh()

    wells = Wells(M, load=load)

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

    return M, elements_lv0, data_impress, wells
