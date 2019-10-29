from ..running.run_simulation import RunSimulation
from .. import directories as direc

state = int(direc.data_loaded['state'])
rodar = RunSimulation(state=state)