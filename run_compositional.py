from packs.directories import data_loaded
from packs import directories as direc
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.compositionalIMPEC import CompositionalFVM
from packs.compositional.stability_check import StabilityCheck
from packs.compositional.properties_calculation import PropertiesCalc
from packs.compositional.update_time import delta_time
from get_inputs_compositional import FluidProperties
from packs.utils import constants as ctes
import os
import numpy as np
import time

class run_simulation:
    '''Class created to compute simulation properties at each simulation time'''
    def __init__(self, name_current, name_all):
        self.name_current_results =os.path.join(direc.flying, name_current + '.npy')
        self.name_all_results = os.path.join(direc.flying, name_all)
        self.loop = 0
        self.vpi = 0.0
        self.t = 0.0
        self.oil_production = 0.
        self.gas_production = 0.
        self.use_vpi = data_loaded['use_vpi']
        self.vpi_save = data_loaded['compositional_data']['vpis_para_gravar_vtk']
        self.time_save = np.array(data_loaded['compositional_data']['time_to_save'])
        self.delta_t = data_loaded['compositional_data']['time_data']['delta_t_ini']
        self.mesh_name =  'compositional_'
        self.all_results = self.get_empty_current_compositional_results()
        self.p1 = PropertiesCalc()

    def initialize(self, load, convert, mesh):
        ''' Function to initialize mesh (preprocess) get and compute initial mesh \
        properties '''
        M, elements_lv0, data_impress, wells = initial_mesh(mesh, load=load, convert=convert)
        ctes.init(M, wells)
        ctes.component_properties()
        fprop = self.get_initial_properties(M, wells)
        return M, data_impress, wells, fprop, load

    def get_initial_properties(self, M, wells):
        ''' get initial fluid - oil, gas and water data and calculate initial \
        properties'''

        fprop = FluidProperties(M, wells) # load reservoir properties data and initialize other data

        '------------------------- Perform initial flash ----------------------'

        if ctes.load_k:
            self.p2 = StabilityCheck(fprop.P, fprop.T)
            fprop.L, fprop.V, fprop.xkj[0:ctes.Nc, 0, :], \
            fprop.xkj[0:ctes.Nc, 1, :], fprop.Csi_j[:,0,:], \
            fprop.Csi_j[:,1,:], fprop.rho_j[:,0,:], fprop.rho_j[:,1,:]  =  \
            self.p2.run_init(fprop.P, fprop.z)

        else: fprop.x = []; fprop.y = []
        if ctes.load_w: fprop.inputs_water_properties(M) #load water properties

        '----------------------- Calculate fluid properties -------------------'

        self.p1.run_outside_loop(M, fprop)

        return fprop

    def run(self, M, wells, fprop, load):
        ''' Function created to compute reservoir and fluid properties at each \
        time step '''

        t0 = time.time()
        t_obj = delta_time(fprop) #get wanted properties in t=n

        '---- Get pressure field and new time step (if the past time step does \
        not obey the CFL condition) -------------------------------------------'

        self.delta_t = CompositionalFVM()(M, wells, fprop, self.delta_t)

        self.t += self.delta_t

        '----------------- Perform Phase stability test and flash -------------'

        if ctes.load_k and ctes.compressible_k:
            fprop.L, fprop.V, fprop.xkj[0:ctes.Nc, 0, :], \
            fprop.xkj[0:ctes.Nc, 1, :], fprop.Csi_j[:,0,:], \
            fprop.Csi_j[:,1,:], fprop.rho_j[:,0,:], fprop.rho_j[:,1,:]  =  \
            self.p2.run(wells, fprop.P, fprop.z)

        '----------------------- Update fluid properties ----------------------'

        self.p1.run_inside_loop(M, fprop)

        '-------------------- Advance in time and save results ----------------'

        self.update_vpi(fprop, wells)
        self.delta_t = t_obj.update_delta_t(self.delta_t, fprop, ctes.load_k, self.loop)#get delta_t with properties in t=n and t=n+1
        if len(wells['ws_p'])>0:self.update_production(fprop, wells)
        self.update_loop()
        t1 = time.time()
        dt = t1 - t0
        if self.use_vpi:
            if np.round(self.vpi,3) in self.vpi_save:
                self.update_current_compositional_results(M, wells, fprop, dt) #ver quem vou salvar
        else:
            if self.time_save[0] == 0.0 or self.t in self.time_save:
                self.update_current_compositional_results(M, wells, fprop, dt)

    def update_loop(self):
        ''' Function to count how many loops it has been since the simulation \
        started'''

        self.loop += 1

    def update_vpi(self, fprop, wells):
        ''' Function to update time in vpi units (volume poroso injetado)'''

        if len(wells['ws_inj'])>0:
            flux_vols_total = wells['values_q_vol']
            flux_total_inj = np.absolute(flux_vols_total)
        else: flux_total_inj = np.zeros(2)

        self.vpi = self.vpi + (flux_total_inj.sum())/sum(fprop.Vp)*self.delta_t

    def get_empty_current_compositional_results(self):
        return [np.array(['loop', 'vpi [s]', 'simulation_time [s]', 't [s]', 'pressure [Pa]', 'Sw', 'So', 'Sg',
                        'Oil_p', 'Gas_p', 'z', 'centroids'])]

    def update_production(self, fprop, wells):
        ''' Function to compute oil and gas production rate [m³/s] through time'''

        self.oil_production +=  abs(fprop.q_phase[:,0,:].sum()) *self.delta_t
        self.gas_production +=  abs(fprop.q_phase[:,1,:].sum())*self.delta_t

    def update_current_compositional_results(self, M, wells, fprop, simulation_time: float = 0.0):

        #total_flux_internal_faces = fprop.total_flux_internal_faces.ravel() #* M.faces.normal[M.faces.internal]
        #total_flux_internal_faces_vector = fprop.total_flux_internal_faces.T * np.abs(M.faces.normal[M.faces.internal])

        self.current_compositional_results = np.array([self.loop, self.vpi, simulation_time,
        self.t, fprop.P, fprop.Sw, fprop.So, fprop.Sg, self.oil_production,
        self.gas_production, fprop.z, M.data['centroid_volumes']],dtype=object)
        self.all_results.append(self.current_compositional_results)
        M.data['saturation'] = fprop.Sw
        M.data['So'] = fprop.So
        M.data['Sg'] = fprop.Sg
        M.data['P'] = fprop.P
        M.data.update_variables_to_mesh()
        M.core.print(file = self.name_all_results + str(self.loop), extension ='.vtk')

    def export_current_compositional_results(self):
         np.save(self.name_current_results, self.current_compositional_results)

    def export_all_results(self):
         np.save(self.name_all_results + str(self.loop) + '.npy', np.array(self.all_results))
         self.all_results = self.get_empty_current_compositional_results()

    def save_infos(self, data_impress, M):
         self.export_current_compositional_results()
         self.export_all_results()
         #data_impress.update_variables_to_mesh()
         #data_impress.export_all_datas_to_npz()
         #M.core.print(file=self.mesh_name, extension='.h5m', config_input="input_cards/print_settings.yml")
        # M.core.print(file = self.mesh_name, extension='.vtk', config_input="input_cards/print_settings.yml")
