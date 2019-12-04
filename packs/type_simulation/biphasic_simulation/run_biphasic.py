from ...simulations.init_simulation import rodar
M = rodar.M
from .biphasic_tpfa import biphasicTpfa
from .biphasic_tpfa import direc

load = direc.data_loaded['load_biphasic']

b1 = biphasicTpfa(M, load=load)

def running():
    verif = True

    while verif:
        b1.run()
