
from ..data_class.data_impress import Data
from ..data_class.elements_lv0 import ElementsLv0
from ..contours.wells import Wells
from ..convert_unit.conversion import Conversion
from ..preprocess.preprocess1 import set_saturation_regions
from ..preprocess.prep0_0 import Preprocess0
from ..directories import data_loaded
from ..multiscale.preprocess.dual_primal.create_dual_and_primal_mesh import MultilevelData
from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
import pdb
import numpy as np

def initial_mesh (mesh, load=False, convert=False):

    compositional = data_loaded['compositional_data']
    load_compositional_data = data_loaded['load_data']

    if compositional and not load_compositional_data:
        M = msh(mesh, dim = 3)
        elements_lv0 = ElementsLv0(M, load=load)
        data_impress = Data(M, elements_lv0, load=load)
        if not load:
            Preprocess0(M, elements_lv0)
        wells = Wells(M, load=load)
        set_saturation_regions(M, wells)
        #fluid_properties = FluidProperties()

    return M,elements_lv0, data_impress, wells #fluid_properties
