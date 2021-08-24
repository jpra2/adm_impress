from packs.data_class import compositional_data


import copy
def update_variables_for_initial_run(fprop, sim, compositional_data2):
    compositional_data = copy.deepcopy(compositional_data2)
    fprop.P = compositional_data['pressure']
    fprop.z = compositional_data['global_composition']
    # fprop.Sg = compositional_data['Sg']
    # fprop.Sw = compositional_data['Sw']
    # fprop.So = compositional_data['So']
    fprop.Nk = compositional_data['mols']
    # fprop.xkj = compositional_data['xkj']
    
    loop_array = compositional_data['loop_array']
    sim.loop = loop_array['loop'][0]
    sim.t = loop_array['t'][0]
    sim.vpi = loop_array['vpi'][0]
    sim.oil_production = loop_array['oil_production'][0]
    sim.gas_production = loop_array['gas_production'][0]
    