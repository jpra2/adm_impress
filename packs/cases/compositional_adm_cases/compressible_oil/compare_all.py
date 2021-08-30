import numpy as np
from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
from packs.cases.compositional_adm_cases.compressible_oil import all_functions
import matplotlib.pyplot as plt

cases_str = ['case3_finescale_3k', 'case9_adm_3k', 'case10_adm_3k', 'case11_adm_3k']
x_legend = ['Finescale', 'Cr=5', 'Cr=10', 'Cr=25']

datas_case = all_functions.load_cases_from_keyword(cases_str, ['loop_array'])

structured_loop_arrays = all_functions.get_loop_array_structured_from_data_cases(datas_case)
sorted_loop_arrays = all_functions.reordenate_loop_arrays_by_loop(structured_loop_arrays)

n_loops = []
total_time_simulation = []
for loop_array in sorted_loop_arrays:
    n_loops.append(loop_array['loop'].shape[0])
    total_time_simulation.append(loop_array['simulation_time'].sum())


fig, ax1 = plt.subplots(1, 1)

# ax1.plot(case1_time, case1_oil_production, '-', label='Finescale')
ax1.bar(x_legend, n_loops, width=0.5, label='N_loops')
ax1.bar(x_legend, total_time_simulation, width=0.3, label='Total_time')
# ax1.set_ylabel('Oil production')
# ax1.set_xlabel('time [days]')
ax1.legend()

plt.savefig('fig_all_comparations' + '.png')
import pdb; pdb.set_trace()