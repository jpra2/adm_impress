import numpy as np
from packs.directories import data_loaded
from packs import directories as direc
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.compositionalIMPEC import CompositionalFVM
from packs.compositional.stability_check import StabilityCheck
from packs.compositional.properties_calculation import PropertiesCalc
from packs.compositional.update_time import delta_time
from update_inputs_compositional import ComponentProperties, FluidProperties
import update_inputs_compositional
from packs.utils import constants as ctes
import os
import time

def initialize(load, convert, mesh):
    M, elements_lv0, data_impress, wells = initial_mesh(mesh, load=load, convert=convert)
    ctes.init(M)
    fprop, kprop = get_initial_properties(M, wells, load, data_loaded)
    return M, data_impress, wells, fprop, kprop, load

def get_initial_properties(M, wells, load, data_loaded):
    kprop = ComponentProperties(data_loaded)
    P = np.array(data_loaded['Pressure']['r1']['value']).astype(float)
    fprop = FluidProperties(kprop)

    if kprop.load_k:
        fprop_block = StabilityCheck(P, ctes.T, kprop)
        fprop_block.run(fprop.z, kprop)
        fprop.run_inputs_k(fprop_block, kprop, ctes.n_volumes)
    else: fprop.x = []; fprop.y = []

    if kprop.load_w: fprop.run_inputs_w(ctes.T, P, data_loaded, ctes.n_volumes)

    PropertiesCalc().run_outside_loop(M, fprop, kprop)
    return fprop, kprop

class run_simulation:

    def __init__(self, name_current, name_all):
        self.name_current_compositional_results =os.path.join(direc.flying, name_current + '.npy')
        self.name_all_compositional_results = os.path.join(direc.flying, name_all)
        self.loop = 0
        self.vpi = 0.0
        self.t = 0.0
        self.use_vpi = data_loaded['use_vpi']
        self.vpi_save = data_loaded['compositional_data']['vpis_para_gravar_vtk']
        self.time_save = data_loaded['compositional_data']['time_to_save']
        self.delta_t = data_loaded['compositional_data']['time_data']['delta_t_ini']
        self.mesh_name =  'compositional_'
        self.all_compositional_results = self.get_empty_current_compositional_results()

    def run(self, M, wells, fprop, kprop, load):
        t0 = time.time()
        t_obj = delta_time(fprop) #get wanted properties in t=n
        self.delta_t = CompositionalFVM().runIMPEC(M, wells, fprop, kprop, self.delta_t)

        self.t += self.delta_t

        if kprop.load_k and kprop.compressible_k:
            for i in range(0, ctes.n_volumes):
                P = fprop.P[i]
                z = fprop.z[0:kprop.Nc,i] #água não entra
                fprop_block = StabilityCheck(P, fprop.T, kprop)
                fprop_block.run(z, kprop)
                fprop.update_all_volumes(fprop_block, i)

        PropertiesCalc().run_inside_loop(M, fprop, kprop)
        self.update_vpi(kprop, fprop, wells)
        self.delta_t = t_obj.update_delta_t(self.delta_t, fprop, kprop.load_k, self.loop)#get delta_t with properties in t=n and t=n+1
        self.update_loop()
        t1 = time.time()
        dt = t1 - t0

        # Talvez isso esteja antes de self.all_compositional_results dentro de update_current_compositional_results
        if self.use_vpi:
            if np.round(self.vpi,3) in self.vpi_save:
                self.update_current_compositional_results(M, wells, fprop, dt) #ver quem vou salvar
        else:
            if np.round(self.t) in self.time_save:
                self.update_current_compositional_results(M, wells, fprop, dt)

    def update_loop(self):
        self.loop += 1

    def update_vpi(self, kprop, fprop, wells):
        if len(wells['ws_inj'])>0:
            flux_vols_total = wells['values_q']
            flux_total_inj = np.absolute(flux_vols_total)
        else: flux_total_inj = np.zeros(2)

        self.vpi = self.vpi + (flux_total_inj.sum())/sum(fprop.Vp)*self.delta_t


    def get_empty_current_compositional_results(self):

        return [np.array(['loop', 'vpi [s]', 'simulation_time [s]', 't [s]', 'pressure [Pa]', 'Sw', 'centroids'])]

    def update_current_compositional_results(self, M, wells, fprop, simulation_time: float = 0.0):

        #total_flux_internal_faces = fprop.total_flux_internal_faces.ravel() #* M.faces.normal[M.faces.internal]
        #total_flux_internal_faces_vector = fprop.total_flux_internal_faces.T * np.abs(M.faces.normal[M.faces.internal])

        self.current_compositional_results = np.array([self.loop, self.vpi, simulation_time,
        self.t, fprop.P, fprop.Sw, M.data['centroid_volumes']])
        self.all_compositional_results.append(self.current_compositional_results)

    def export_current_compositional_results(self):
         np.save(self.name_current_compositional_results, self.current_compositional_results)

    def export_all_compositional_results(self):
         np.save(self.name_all_compositional_results + str(self.loop) + '.npy', np.array(self.all_compositional_results))
         self.all_compositional_results = self.get_empty_current_compositional_results()

    def save_infos(self, data_impress, M):
         self.export_current_compositional_results()
         self.export_all_compositional_results()
         data_impress.update_variables_to_mesh()
         data_impress.export_all_datas_to_npz()
         M.core.print(file=self.mesh_name, extension='.h5m', config_input="input_cards/print_settings.yml")
         # M.core.print(file=self.mesh_name, extension='.vtk', config_input="input_cards/print_settings.yml")
