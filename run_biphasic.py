from packs.simulations.init_simulation import rodar
from packs.type_simulation.biphasic_simulation.biphasic_tpfa import biphasicTpfa
from packs import directories as direc
import os
import numpy as np


load = np.load(direc.name_load)[0]
verif = True
M = rodar.M
b1 = biphasicTpfa(M, load=load)

while verif:
    if b1.loop % b1.loops_para_gravar == 0 and b1.loop > 0:
        b1.run(save=True)
        verif = False

    else:
        b1.run()
