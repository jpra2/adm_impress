from packs.biphasic.biphasic_tpfa import BiphasicTpfa
from packs.running.initial_mesh_properties import initial_mesh
from packs.utils.info_manager import InfoManager
import os

dd = InfoManager('input_cards/inputs0_2.yml', 'input_cards/inputs0.yml')
dd.save_obj()
# import pdb; pdb.set_trace()
load = dd['load_data']
convert = dd['convert_english_to_SI']
M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
b1 = BiphasicTpfa(data_impress, elements_lv0, wells)
# b1.run()
n = 5
cont = 1
verif = True
while verif:
    b1.run()
    if cont % n == 0:
        cont = 1
        data_impress.update_variables_to_mesh()
        name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
        M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")
        import pdb; pdb.set_trace()
    cont += 1


# b1.run()
# b1.run()
# b1.run()
# b1.run()
# b1.run()
# b1.run()
# b1.run()
data_impress.update_variables_to_mesh()
name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")


import pdb; pdb.set_trace()
