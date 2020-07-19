from packs.tpfa.biphasic.values_biphasic import BiphasicData
from packs.tpfa.biphasic.get_relative import get_relative_permeability
from packs.tpfa.biphasic.biphasic_fluid_properties import BiphasicFluidProperties
from packs.tpfa.biphasic.biphasic_utils import BiphasicUtils
from packs.tpfa.common_files import SimulationVariables

class TpfaBiphasicCons():

    def __init__(self, load=False):
        self.data = BiphasicData(load=load)
        self.relative_permeability = get_relative_permeability()
        self.properties = BiphasicFluidProperties()
        self.biphasic_utils = BiphasicUtils()
        self.simulation_variables = SimulationVariables()

    def set_current_saturation(self, saturation):
        self.data['saturation'] = saturation

    def get_empty_current_biphasic_results(self):

        dty = [('loop', np.int), ('delta_t [s]', np.float), ('simulation_time [s]', np.float),
               ('oil_production [m3/s]', np.float), ('water_production [m3/s]', np.float),
               ('t [s]', np.float), ('wor', np.float), ('vpi', np.float), ('contador_vtk', np.int)]

        # return [np.array(['loop', 'delta_t [s]', 'simulation_time [s]',
        #     'oil_production [m3/s]', 'water_production [m3/s]', 't [s]', 'wor', 'vpi', 'contador_vtk'])]

        return 0
