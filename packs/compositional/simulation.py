import numpy as np
from ..directories import data_loaded
from .compositionalIMPEC import CompositionalFVM
from .stability_check import StabilityCheck
from .properties_calculation import PropertiesCalc
from .update_time import delta_time
import update_inputs_compositional
import time

class run_simulation:
    def __init__(self):
        self.name_current_compositional_results = 'flying/results_caso.npy'
        self.loop = 0
        self.vpi = 0.0
        self.t = 0.0
        self.contador_vtk = 0
        self.deltaT = 0.002 # chute? talvez nao

    def run(self, M, data_impress, wells, fprop, fprop_block, kprop, load, n_volumes):
        t0 = time.time()
        t_obj = delta_time(fprop) #get wanted properties in t=n
        self.t += self.deltaT

        ss = CompositionalFVM(M, data_impress, wells, fprop, fprop_block, kprop, load, self.deltaT)
        prop = PropertiesCalc(data_impress, wells, fprop)
        prop.run_inside_loop(data_impress, wells, fprop)

        for i in range(1, n_volumes):
            P = fprop.P[i]
            z = fprop.z[0:fprop.Nc,i] #água não entra
            fprop_block = StabilityCheck(z, P, fprop.T, fprop.R, fprop.Nc, kprop)
            fprop_block.run(z, kprop)
            fprop.update_all_volumes(fprop_block, i)

        self.deltaT = t_obj.update_deltaT(data_loaded, wells, self.deltaT, fprop)#get deltaT with properties in t=n and t=n+1
        self.update_loop()
        t1 = time.time()
        dt = t1 - t0
        # self.update_current_compositional_results(wells, dt) #ver quem vou salvar
    def update_loop(self):
        self.loop += 1
    #
    # def get_empty_current_compositional_results(self):
    #
    #     return [np.array(['loop', 'delta_t [s]', 'simulation_time [s]',
    #         'oil_production [m3/s]', 'water_production [m3/s]', 't [s]', 'wor', 'vpi', 'contador_vtk'])]
    #
    # def update_current_compositional_results(self, wells, data_impress, simulation_time: float = 0.0):
    #
    #     water_production = sum(flux_vols_total[fprop.Nc, wells[ws_prod]])
    #     oil_production = (flux_vols_total[0:fprop.Nc, wells[ws_prod]]).sum(axis = 0).sum(axis = 0)
    #
    #     wor = water_production / oil_production
    #
    #     self.current_compositional_results = np.array([self.loop, self.delta_t, simulation_time,
    #         -oil_production, -water_production, self.t, wor, self.vpi, self.contador_vtk])
    #
    #     self.all_compositional_results.append(self.current_compositional_results)
    #
    # def export_current_biphasic_results(self):
    #     np.save(self.name_current_compositional_results, self.current_compositional_results)
    #
    # def export_all_biphasic_results(self):
    #     np.save(self.name_all_compositional_results + str(self.loop) + '.npy', np.array(self.all_compositional_results))
    #     self.all_compositional_results = self.get_empty_current_biphasic_results()
    #
    # def save_infos(self, data_impress):
    #     self.export_current_biphasic_results()
    #     self.export_all_biphasic_results()
    #     data_impress.update_variables_to_mesh()
    #     data_impress.export_all_datas_to_npz()
    #     self.mesh.core.print(file=self.mesh_name, extension='.h5m', config_input="input_cards/print_settings.yml")
