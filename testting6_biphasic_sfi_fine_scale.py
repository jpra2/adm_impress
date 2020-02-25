from packs.biphasic.biphasic_sfi import BiphasicSfi
from packs.running.initial_mesh_properties import initial_mesh
import os

M, elements_lv0, data_impress, wells = initial_mesh()
b1 = BiphasicSfi(M, data_impress, elements_lv0, wells)

n = 10
n2 = 20
n_for_save = 20
cont_for_save = 1
loop = 0
cont = 1
cont2 = 1
verif = True

# import pdb; pdb.set_trace()

while verif:
    if cont_for_save % n_for_save == 0:
        save = True
        cont_for_save = 1
    else:
        save = False
        cont_for_save += 1
    resp = 1
    while resp != 0:
        b1.run()
        resp = b1.run_2(save=save)

    import pdb; pdb.set_trace()
    print(f'\n loop: {b1.loop}\n')
    if cont % n == 0:
        cont = 1
        data_impress.update_variables_to_mesh()
        name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
        M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")
        # import pdb; pdb.set_trace()
    cont += 1
    loop += 1
    if loop > 200:
        if cont2 % n2 == 0:
            cont2 = 1
            data_impress.update_variables_to_mesh()
            name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
            M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")
            # import pdb; pdb.set_trace()

        cont2 += 1

import pdb; pdb.set_trace()
