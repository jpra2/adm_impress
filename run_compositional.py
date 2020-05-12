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
import os
import time

def initialize(load, convert, mesh):
    M, elements_lv0, data_impress, wells = initial_mesh(mesh, load=load, convert=convert)
    n_volumes = data_impress.len_entities['volumes']
    fprop, kprop = get_initial_properties(M, data_impress, wells, load, data_loaded, n_volumes)
    return M, data_impress, wells, fprop, kprop, load, n_volumes

def get_initial_properties(M, data_impress, wells, load, data_loaded, n_volumes):
    kprop = ComponentProperties(data_loaded)
    P, T = update_inputs_compositional.inputs_overall_properties(data_loaded)
    fprop = FluidProperties(kprop)
    get_constant_properties(data_impress, fprop)

    fprop_block = StabilityCheck(P, T, kprop)
    if kprop.load_k:
        fprop_block.run(fprop.z, kprop)
        fprop.run_inputs_k(fprop_block, kprop, n_volumes)
    else: fprop.x = []; fprop.y = []
    if kprop.load_w: fprop.run_inputs_w(T, P, data_loaded, n_volumes)

    PropertiesCalc(n_volumes).run_outside_loop(data_impress, wells, fprop, kprop)
    return fprop, kprop

def get_constant_properties(data_impress, fprop):
    fprop.Pf = np.array(data_loaded['compositional_data']['Pf']).astype(float)
    fprop.porosity = data_impress['poro']
    fprop.cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)
    fprop.Pw = np.array(data_loaded['compositional_data']['water_data']['Pw']).astype(float)
    fprop.Cw = np.array(data_loaded['compositional_data']['water_data']['Cw']).astype(float)
    fprop.Vbulk = data_impress['volume']

class run_simulation:

    def __init__(self, M, wells, data_impress, name_current, name_all):
        self.name_current_compositional_results =os.path.join(direc.flying, name_current + '.npy')
        self.name_all_compositional_results = os.path.join(direc.flying, name_all)
        self.loop = 0
        self.vpi = 0.0
        self.t = 0.0
        self.contador_vtk = 0
        self.n_volumes = len(data_impress['volume'])
        self.use_vpi = data_loaded['use_vpi']
        self.vpi_save = data_loaded['compositional_data']['vpis_para_gravar_vtk']
        self.time_save = data_loaded['compositional_data']['time_to_save']
        self.delta_t = data_loaded['compositional_data']['time_data']['delta_t_ini']
        self.mesh_name =  'compositional_'
        self.all_compositional_results = self.get_empty_current_compositional_results()
        self.FVM = CompositionalFVM(M, self.n_volumes, wells)

    def run(self, M, data_impress, wells, fprop, kprop, load, n_volumes):
        t0 = time.time()
        t_obj = delta_time(fprop) #get wanted properties in t=n
        #import pdb; pdb.set_trace()
        self.delta_t = self.FVM.runIMPEC(M, data_loaded, data_impress, wells, fprop, kprop, self.delta_t)

        self.t += self.delta_t

        if kprop.load_k and kprop.compressible_k:
            for i in range(0, n_volumes):
                P = fprop.P[i]
                z = fprop.z[0:kprop.Nc,i] #água não entra
                fprop_block = StabilityCheck(P, fprop.T, kprop)
                fprop_block.run(z, kprop)
                fprop.update_all_volumes(fprop_block, i)

        PropertiesCalc(n_volumes).run_inside_loop(data_impress, wells, fprop, kprop)
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
        #self.phase_existance = np.zeros([1, kprop.n_phases, self.n_volumes])

        #if kprop.load_k:
        #    self.phase_existance[0,0,:] = np.sign(fprop.L)
        #    self.phase_existance[0,1,:] = np.sign(fprop.V)**2
        #if kprop.load_w:
        #    self.phase_existance[0,kprop.n_phases-1,:] = 1

        if len(wells['ws_inj'])>0:
            #flux_vols_total = fprop.component_flux_vols_total / \
            #(fprop.component_molar_fractions * self.phase_existance * fprop.phase_molar_densities).sum(axis=1)[:,wells['ws_inj']]

            flux_vols_total = wells['values_q']
            # (fprop.component_molar_fractions * phase_existance * fprop.phase_molar_densities).sum(axis=1)[:,wells['ws_inj']]
            flux_total_inj = np.absolute(flux_vols_total)
        else: flux_total_inj = np.zeros(2)
        self.vpi += (flux_total_inj.sum()*self.delta_t)/sum(fprop.Vp)


    def get_empty_current_compositional_results(self):

        return [np.array(['loop', 'delta_t [s]', 'simulation_time [s]', 't [s]', 'pressure [Pa]', 'Sw', 'centroids'])]

    def update_current_compositional_results(self, M, wells, fprop, simulation_time: float = 0.0):

        #total_flux_internal_faces = fprop.total_flux_internal_faces.ravel() #* M.faces.normal[M.faces.internal]
        #total_flux_internal_faces_vector = fprop.total_flux_internal_faces.T * np.abs(M.faces.normal[M.faces.internal])

        self.current_compositional_results = np.array([self.loop, self.delta_t, simulation_time,
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
