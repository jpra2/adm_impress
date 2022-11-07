import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy import interpolate


flying = 'flying'
name = 'results_'
arquivos = os.listdir(flying)

for  arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_Hoteit_Firoo_2k_ex1_vpi_IMPEC_FOU_CFL_0_9_3936.npy', allow_pickle=True)

        for data in datas[1:]:

            import pdb; pdb.set_trace()
