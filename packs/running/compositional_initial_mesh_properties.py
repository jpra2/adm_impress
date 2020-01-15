
from ..data_class.data_impress import Data
from ..data_class.elements_lv0 import ElementsLv0
from ..contours.wells import Wells
from ..convert_unit.conversion import Conversion
from ..preprocess.preprocess1 import set_saturation_regions
from ..preprocess.prep0_0 import Preprocess0
from ..directories import data_loaded
from ..multiscale.preprocess.dual_primal.create_dual_and_primal_mesh import MultilevelData
import numpy as np

def initial_mesh(load=False, convert=False):

    global data_loaded

    compositional = data_loaded['compositional_data']
    load_compositional_data = data_loaded['load_data']

    if compositional and not load_compositional_data:
        #set_saturation_regions(M, wells)
        # TODO: atualizar gama

    return M,elements_lv0,data_impress,wells,fluid_properties
