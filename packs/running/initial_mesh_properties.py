
from ..data_class.data_impress import Data
from ..data_class.elements_lv0 import ElementsLv0
from ..contours.wells import Wells
from ..convert_unit.conversion import Conversion
from ..preprocess.preprocess1 import set_saturation_regions
from ..preprocess.prep0_0 import Preprocess0
from ..directories import data_loaded

# import pdb; pdb.set_trace()

def initial_mesh(load=False, convert=False):
    from ..load.preprocessor0 import M
    global data_loaded

    elements_lv0 = ElementsLv0(M, load=load)
    data_impress = Data(M, elements_lv0, load=load)
    if not load:
        Preprocess0(M, elements_lv0)

    wells = Wells(M, load=load)

    biphasic = data_loaded['biphasic']
    if biphasic:
        set_saturation_regions(M, wells)

    if convert:
        conversion = Conversion(wells, data_impress)
        conversion.convert_English_to_SI()

    if not load:

        wells.update_values_to_mesh()
        data_impress.update_variables_to_mesh()
        wells.export_all_datas_to_npz()

    return M, elements_lv0, data_impress, wells
