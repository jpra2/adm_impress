import copy
import numpy as np
import yaml

class TimeProfile:
    dt_update_boundary_nodes_weights_neumann = 0
    dt_update_boundary_nodes_weights_neumann_mean = 0
    dt_set_local_neumann_problem = 0.0
    dt_weight_finescale = 0.0
    dt_update_coarse_weight = 0.0
    dt_update_coarse_weight_mean = 0.0
    percent_weight_update_in_neumann = 0.0
    dt_solution_ms = 0.0
    dt_neumann = 0.0
    n_coarses = 0.0
    n_pressure_updates = 0.0
    dt_total_pressure_update = 0.0
    percent_neumann_in_total_pressure = 0.0
    dt_set_adm_mesh = 0.0
    dt_set_finescale_problem = 0.0
    dt_solution_fs = 0.0
    dt_update_mobility = 0.0
    dt_local_flux_update = 0.0
    dt_neumann_preprocess = 0.0
    
    
    @classmethod
    def update_data(cls):
        ## cls.dt_update_coarse_weight_mean eh uma parte de dt_update_boundary_nodes_weights_neumann_mean
        cls.dt_update_boundary_nodes_weights_neumann_mean = cls.dt_update_boundary_nodes_weights_neumann/cls.n_coarses
        cls.dt_update_coarse_weight_mean = cls.dt_update_coarse_weight/cls.n_coarses
        cls.percent_weight_update_in_neumann = cls.dt_update_boundary_nodes_weights_neumann/cls.dt_neumann
        cls.n_pressure_updates += 1
        cls.percent_neumann_in_total_pressure = cls.dt_neumann/cls.dt_total_pressure_update
    
    @classmethod
    def export_data_dict(cls):
        
        exclude_names = np.array(['__module__', 'update_data', 'export_data_dict', '__dict__',
                         '__weakref__', '__doc__', 'show_data', 'reset_data', 'export_to_data'])
        all_names = np.array(list(cls.__dict__.keys()))
        my_names = np.setdiff1d(all_names, exclude_names)
        
        data = dict()
        for key in my_names:
            data.update({key: cls.__dict__[key]})
        
        return data
    
    @classmethod    
    def show_data(cls):
        data = cls.export_data_dict()
        for key in data.keys():
            print(f"{key}: {data[key]}")
    
    @classmethod
    def reset_data(cls):
        cls.dt_update_boundary_nodes_weights_neumann = 0.0
        cls.dt_set_local_neumann_problem = 0.0
        cls.dt_local_flux_update = 0.0
        cls.dt_neumann_preprocess = 0.0
    
    @classmethod
    def export_to_data(cls, filename_export_times:str, **kwargs):
        data = cls.export_data_dict()
        with open(filename, 'w') as f:
            yaml.dump(data, f)
        
        
        
        
        
        
        
        
        
    
    