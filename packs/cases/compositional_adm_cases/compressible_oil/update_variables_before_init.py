import copy
from packs.data_class.sparse_operators import SparseOperators
from packs.data_class.compositional_data import CompositionalData

def get_copy(data):
    return copy.deepcopy(data)

def update_variables_for_initial_run_adm(fprop, sim, latest_mobilities, compositional_data2, OP_AMS, manage_operators2):
    compositional_data: CompositionalData = compositional_data2
    manage_operators: SparseOperators = manage_operators2
    compositional_data.load_from_npz()
    manage_operators.load()
    fprop.P = get_copy(compositional_data['pressure'])
    fprop.z = get_copy(compositional_data['global_composition'])
    # fprop.Sg = compositional_data['Sg']
    # fprop.Sw = compositional_data['Sw']
    # fprop.So = compositional_data['So']
    fprop.Nk = get_copy(compositional_data['mols'])
    fprop.Vp = get_copy(compositional_data['Vp'])
    # fprop.xkj = compositional_data['xkj']
    latest_mobilities[:] = get_copy(compositional_data2['latest_mobility'])
    
    loop_array = compositional_data['loop_array']
    sim.loop = loop_array['loop'][0]
    sim.t = loop_array['t'][0]
    sim.vpi = loop_array['vpi'][0]
    sim.oil_production = loop_array['oil_production'][0]
    sim.gas_production = loop_array['gas_production'][0]
    
    OP_AMS[:] = manage_operators['prolongation_level_1']

def update_variables_for_initial_run_finescale(fprop, sim, compositional_data2):
    compositional_data: CompositionalData = compositional_data2
    compositional_data.load_from_npz()
    fprop.P = get_copy(compositional_data['pressure'])
    fprop.z = get_copy(compositional_data['global_composition'])
    fprop.q = get_copy(compositional_data['q'])
    # fprop.Sg = compositional_data['Sg']
    # fprop.Sw = compositional_data['Sw']
    # fprop.So = compositional_data['So']
    fprop.Nk = get_copy(compositional_data['mols'])
    fprop.Vp = get_copy(compositional_data['Vp'])
    # fprop.xkj = compositional_data['xkj']
    
    loop_array = compositional_data['loop_array']
    sim.loop = loop_array['loop'][0]
    sim.t = loop_array['t'][0]
    sim.vpi = loop_array['vpi'][0]
    sim.oil_production = loop_array['oil_production'][0]
    sim.gas_production = loop_array['gas_production'][0]