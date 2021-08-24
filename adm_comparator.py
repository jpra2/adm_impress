from run_test_cases_mono import *
import time

t0=time.time()

remove_previous_files()
run_test_cases()
all_cases_results=organize_results()

# print_results_single(all_cases_results)
print_results_two_phase(all_cases_results)
print(time.time()-t0, 'total full time!!!')
