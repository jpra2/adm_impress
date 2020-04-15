from packs.biphasic.biphasic_tpfa import BiphasicTpfa
from packs.running.initial_mesh_properties import initial_mesh
from packs.utils.info_manager import InfoManager
from packs.directories import data_loaded as dd
import os

M, elements_lv0, data_impress, wells = initial_mesh()
b1 = BiphasicTpfa(M, data_impress, elements_lv0, wells)
# b1.run()
# import pdb; pdb.set_trace()
n = 1
n2 = 20
n_for_save = 20
cont_for_save = 1
loop = 0
cont = 1
cont2 = 1
verif = True

meshset_volumes = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
# import pdb; pdb.set_trace()
while verif:
    if cont_for_save % n_for_save == 0:
        b1.run(save=True)
        cont_for_save = 1
    else:
        b1.run()
        cont_for_save += 1
    print(f'\n loop: {b1.loop}\n')
    if cont % n == 0:
        cont = 1
        data_impress.update_variables_to_mesh()
        name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
        M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")
        # import pdb; pdb.set_trace()
    M.core.mb.write_file('results/testt_'+str(cont)+'.vtk', [meshset_volumes])
    cont += 1
    if cont % n2 == 0:
        import pdb; pdb.set_trace()
    loop += 1
    if loop > 200:
        if cont2 % n2 == 0:
            cont2 = 1
            data_impress.update_variables_to_mesh()
            name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
            M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")
            # import pdb; pdb.set_trace()

        cont2 += 1



# b1.run()
# b1.run()
# b1.run()
# b1.run()
# b1.run()
# b1.run()
# b1.run()
# data_impress.update_variables_to_mesh()
# name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
# M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")


import pdb; pdb.set_trace()
