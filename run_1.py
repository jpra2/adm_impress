
import os
import numpy as np
os.system('python3 update_inputs.py')
# os.system('python3 testting2_biphasic.py')
# os.system('python3 testting1_monophasic_multilevel.py')
os.system('python3 testting3_adm_method.py')
# os.system('python3 testting1_monophasic.py')

# from packs.load.preprocessor0 import M
# # M.pressure[:] = 1.0
# M.pressure[[1,2,3]] = [2.0, 1.5, 4.0]
# M.pressure.update()
# M.permeability[[1,2]] = np.array([
# np.array([1.0, 0.0, 2.0, 0.0, 4.0, 0.0, 8.0, 0.0, 16.0]),
# np.array([5.0, 0.0, 2.0, 0.0, 4.0, 0.0, 8.0, 0.0, 16.0])
# ])
# M.permeability.update()
#
# import pdb; pdb.set_trace()
