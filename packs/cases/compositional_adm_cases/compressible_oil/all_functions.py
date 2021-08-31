from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
import numpy as np

def organize_cases_by_loop(datas, loops):
    loops2 = pd.Series(loops)
    loops2_sorted = loops2.sort_values()
    index_sorted = loops2_sorted.index.values
    resp = np.array(datas)[index_sorted]   
    return resp

def extrair_dado(data_case, key):
    resp = []
    for data in data_case:
        resp.append(data[key])
    
    resp = np.array(resp)
    return resp

def erro_abs(v1, v2):
    return np.absolute(v1 - v2)

def get_data_from_loop_array(keyword, data_case):
    data = []
    for sim in data_case:
        data.append(sim['loop_array'][keyword][0])
    
    return np.array(data)

def create_loop_array_structured(data_case):
    resp = []
    for data in data_case:
        resp.append(data['loop_array'])
    
    return np.array(resp)

def load_cases_from_keyword(cases_str, keywords):
    all_cases = []
    for case_str in cases_str:
        case = CumulativeCompositionalDataManager(description=case_str)
        data_case = case.load_all_datas_from_keys(keywords)
        all_cases.append(data_case)
    
    return all_cases

def get_loop_array_structured_from_data_cases(data_cases):
    loops_array_structured = []
    for data_case in data_cases:
        loops_array_structured.append(create_loop_array_structured(data_case))
    
    return loops_array_structured

def reordenate_loop_arrays_by_loop(loop_arrays):
    resp = []
    for loop_array in loop_arrays:
        resp.append(np.sort(loop_array, order='loop'))
    
    return resp

def get_empty_loop_array():
    loop_array = np.zeros(
        1,
        dtype=[
            ('loop', int),
            ('t', float),
            ('vpi', float),
            ('simulation_time', float),
            ('oil_production', float),
            ('gas_production', float),
            ('n_volumes_update_base_functions', int),
            ('total_volumes_updated', int),
            ('active_volumes', int),
            ('oil_rate', float),
            ('gas_rate', float),
            ('total_simulation_time', float),
            ('n_total_loops', int)   
        ]
    )
    return loop_array
    