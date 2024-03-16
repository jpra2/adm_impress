import numpy as np
from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
from packs.cases.compositional_adm_cases.compressible_oil.all_functions import organize_cases_by_loop, extrair_dado, erro_abs, get_data_from_loop_array, create_loop_array_structured
from packs.cases.compositional_adm_cases.compressible_oil import descriptions

MY_LOOP = 4000
key_list = ['loop_array']
# M = msh('mesh/80x80_BL.msh')
description2 = descriptions.case2_adm_description

case2 = CumulativeCompositionalDataManager(description=description2)
datas_case2 = case2.load_all_datas_from_keys(key_list)

loops2 = []
for data in datas_case2:
    loops2.append(data['loop_array']['loop'][0])
loops2 = np.array(loops2)

datas_case2 = organize_cases_by_loop(datas_case2, loops2)
case2_structured_loop_array = create_loop_array_structured(datas_case2)

case2_time = case2_structured_loop_array['t'].flatten()
case2_time = case2_time/86400

datas = case2_structured_loop_array[case2_structured_loop_array['loop'] == MY_LOOP]

time_days = datas['t']/86400
print(time_days)
print(datas['vpi'])




import pdb; pdb.set_trace()