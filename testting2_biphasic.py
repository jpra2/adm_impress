from packs.biphasic.biphasic_tpfa import BiphasicTpfa
from packs.running.initial_mesh_properties import initial_mesh
from packs.utils.info_manager import InfoManager
from packs.directories import data_loaded as dd
import os
import time
import numpy as np
M, elements_lv0, data_impress, wells = initial_mesh()
b1 = BiphasicTpfa(M, data_impress, elements_lv0, wells)
def export_finescale_results(vals_vpi,vals_delta_t,vals_wor, t_comp):
    vals_vpi=np.array(vals_vpi)
    vals_delta_t=np.array(vals_delta_t)
    vals_wor=np.array(vals_wor)
    t_comp=np.array(t_comp)
    vars=[vpi,vals_delta_t,vals_wor, t_comp]
    names=['vpi','delta_t', 'wor', 't_comp']
    for i in range(len(vars)):
        np.save('results/biphasic/finescale/'+names[i]+'.npy',vars[i])
# b1.run()
# import pdb; pdb.set_trace()
n = 1
n2 = 100
n_for_save = 100
cont_for_save = 1
loop = 0
cont = 1
cont2 = 1
verif = True
vpis_for_save=np.arange(0.0,0.051,0.01)
count_save=0
pp=100
meshset_volumes = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
# import pdb; pdb.set_trace()
vpi=[]
delta_t=[]
wor=[]
t_comp=[]

while verif:
    t00=time.time()
    #
    # if cont_for_save % n_for_save == 0:
    #     b1.run(save=True)
    #     cont_for_save = 1
    # else:
    #     b1.run()
    #     cont_for_save += 1
    b1.run()
    t_comp.append(b1.t_comp)
    vpi.append(b1.vpi)
    delta_t.append(b1.delta_t)
    wor.append(b1.wor)
    print(f'\n loop: {b1.loop}\n')
    # if cont % n == 0:
    #     cont = 1
    #     data_impress.update_variables_to_mesh()
    #     name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
    #     M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")
    #     # import pdb; pdb.set_trace()
    # if cont % pp == 0:
    if vpis_for_save[count_save]<b1.vpi:
        data_impress.update_variables_to_mesh()
        # M.core.mb.write_file('results/testt_'+str(cont)+'.vtk', [meshset_volumes])
        file_count=str(int(100*vpis_for_save[count_save]))

        M.core.mb.write_file('results/biphasic/finescale/vtks/volumes_'+file_count+'.vtk', [meshset_volumes])
        count_save+=1

    if cont % n2 == 0:
        export_finescale_results(vpi, delta_t, wor,t_comp)
        verif=False
    cont += 1
    # loop += 1
    # if loop > 200:
    #     if cont2 % n2 == 0:
    #         cont2 = 1
    #         data_impress.update_variables_to_mesh()
    #         name = os.path.join('results', 'biphasic') + '_loop_' + str(b1.loop)
    #         M.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")
    #         # import pdb; pdb.set_trace()
    #
    #     cont2 += 1
    # print(time.time()-t00,'loop time finescale')
