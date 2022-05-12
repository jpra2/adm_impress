from packs.directories import data_loaded
from packs import directories as direc
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.update_time import delta_time
from get_inputs_compositional import FluidProperties
from packs.utils import constants as ctes
import os
import numpy as np
import time

if data_loaded['compositional_data']['solver']['IMPSAT']:
    from packs.compositional.IMPSAT.compositionalIMPSAT import CompositionalFVM
    from packs.compositional.IMPSAT.properties_calculation import PropertiesCalc
else:
    from packs.compositional.IMPEC.compositionalIMPEC import CompositionalFVM
    from packs.compositional.IMPEC.properties_calculation import PropertiesCalc

if data_loaded['compositional_data']['component_data']['constant_K']:
    from packs.compositional.Kflash import StabilityCheck
else:
    from packs.compositional.stability_check import StabilityCheck

class run_simulation:
    '''Class created to compute simulation properties at each simulation time'''
    def __init__(self, name_current, name_all):
        #name_all = name_all + data_loaded['solver'] + hiperbolic_method + ctes.RS + ctes.
        self.name_current_results =os.path.join(direc.flying, name_current + '.npy')
        self.name_all_results = os.path.join(direc.flying, name_all)
        self.loop = 0
        self.vpi = 0.0
        self.t = 0.0
        self.oil_production = 0.
        self.gas_production = 0.
        self.oil_production_rate = []
        self.gas_production_rate = []
        self.vector_time = []
        self.use_vpi = data_loaded['use_vpi']
        self.vpi_save = data_loaded['compositional_data']['vpis_para_gravar_vtk']
        self.time_save = np.array(data_loaded['compositional_data']['time_to_save'])
        self.delta_t = data_loaded['compositional_data']['time_data']['delta_t_ini']
        self.mesh_name =  'compositional_'
        self.all_results = self.get_empty_current_compositional_results()
        self.p1 = PropertiesCalc()
        self.sim_time = 0

    def initialize(self, load, convert, mesh):
        ''' Function to initialize mesh (preprocess) get and compute initial mesh \
        properties '''
        M, elements_lv0, data_impress, wells = initial_mesh(mesh, load=load, convert=convert)
        ctes.init(M, wells)
        ctes.component_properties()

        if ctes.FR:
            from packs.compositional import prep_FR as ctes_FR
            ctes_FR.run(M)

        if ctes.MUSCL:
            from packs.compositional import prep_MUSCL as ctes_MUSCL
            ctes_MUSCL.run(M)

        fprop = self.get_initial_properties(M, wells)
        return M, data_impress, wells, fprop, load

    def get_well_inj_properties(self, M, fprop, wells):

        z = (wells['z']).T

        #q_wells = wells['ws_inj']#[wells['inj_cond']=='reservoir']
        if ctes.load_k and any(z[0:ctes.Nc].flatten()>0):
            p_well = StabilityCheck(fprop.P[wells['ws_inj']], fprop.T)
            L, V, x, y, Csi_L, Csi_V, rho_L, rho_V  =  \
                p_well.run_init(fprop.P[wells['ws_inj']], z[0:ctes.Nc])

            'if the injector well has a prescribed pressure condition'

            if any(wells['inj_cond']=='reservoir'):
                ws_p_inj_ind = np.argwhere(wells['ws_p']==wells['ws_inj'])
                ws_p_inj = wells['ws_p'][ws_p_inj_ind].flatten()
                ws_inj_p = np.argwhere(wells['ws_inj']==wells['ws_p']).flatten()

                fprop.z[..., ws_p_inj] = z[0:ctes.Nc,ws_inj_p]
                xkj_ws = np.ones((ctes.n_components,ctes.n_phases,len(ws_p_inj)))
                xkj_ws[:ctes.Nc,0] = x[...,ws_inj_p]
                xkj_ws[:ctes.Nc,1] = y[...,ws_inj_p]
                if ctes.load_w:
                    xkj_ws[:ctes.Nc,-1] = 0
                    xkj_ws[-1,:-1] = 0

                mobility_ratio_ws = np.empty((1,ctes.n_phases,len(ws_p_inj)))
                mobility_ratio_ws[:,0] = (L/(L+V))[ws_p_inj]
                mobility_ratio_ws[:,1] = (V/(L+V))[ws_p_inj]
                Csi_j_ws = np.empty((1,ctes.n_phases,len(ws_p_inj)))
                Csi_j_ws[:,0,:] = Csi_L[ws_p_inj]
                Csi_j_ws[:,1,:] = Csi_V[ws_p_inj]
                wells['inj_p_term'] = xkj_ws * mobility_ratio_ws * Csi_j_ws
                'if the injector well has a prescribed flux condition'
                if len(wells['ws_q'])>0:
                    self.q_vol = np.copy(wells['values_q'][:,wells['inj_cond']=='reservoir'])
                    #rever esse Csi_L*L + Csi_V*V
                    wells['values_q'][:,wells['inj_cond']=='reservoir'] = (Csi_V * V + Csi_L * L) * self.q_vol
                    wells['values_q_vol'][:,wells['inj_cond']=='reservoir'] = self.q_vol
            else:
                wells['inj_p_term'] = []
                qk_molar = wells['values_q'][:,wells['inj_cond']=='surface']
                wells['values_q_vol'][:,wells['inj_cond']=='surface'] = qk_molar / \
                    ((Csi_V * V + Csi_L * L))[wells['inj_cond']=='surface']

            ctes.P_SC *= np.ones_like(wells['ws_prod'])
            ctes.T_SC = fprop.T#* np.ones_like(wells['ws_prod'])
            #p_well = StabilityCheck(ctes.P_SC, ctes.T_SC)
            #L, V, x, y, Csi_L, Csi_V, rho_L, rho_V  =  \
            #    p_well.run_init(ctes.P_SC, fprop.z[0:ctes.Nc,wells['ws_prod']])
        else: wells['inj_p_term'] = []

    def get_initial_properties(self, M, wells):
        ''' get initial fluid - oil, gas and water data and calculate initial \
        properties'''

        fprop = FluidProperties(M, wells) # load reservoir properties data and initialize other data

        '------------------------- Perform initial flash ----------------------'
        self.get_well_inj_properties(M, fprop, wells)
        #fprop.z[:,0] = np.array([0,1])
        if ctes.load_k:

            self.p2 = StabilityCheck(fprop.P, fprop.T)
            fprop.L, fprop.V, fprop.xkj[0:ctes.Nc, 0, :], \
            fprop.xkj[0:ctes.Nc, 1, :], fprop.Csi_j[:,0,:], \
            fprop.Csi_j[:,1,:], fprop.rho_j[:,0,:], fprop.rho_j[:,1,:]  =  \
            self.p2.run_init(fprop.P, np.copy(fprop.z))

            ctes.P_SC *= np.ones_like(wells['ws_prod'])
            ctes.T_SC = fprop.T#* np.ones_like(wells['ws_prod'])
            p_well = StabilityCheck(ctes.P_SC, ctes.T_SC)
            L, V, x, y, Csi_L, Csi_V, rho_L, rho_V  =  \
                p_well.run_init(ctes.P_SC, fprop.z[0:ctes.Nc,wells['ws_prod']])

        else: fprop.x = []; fprop.y = []; fprop.L = []; fprop.V = []; wells['inj_term'] = []

        if ctes.load_w:
            fprop.inputs_water_properties(M) #load water properties

        '----------------------- Calculate fluid properties -------------------'

        self.p1.run_outside_loop(M, fprop, wells)
        return fprop

    def run(self, M, wells, fprop, load):
        ''' Function created to compute reservoir and fluid properties at each \
        time step '''

        t0 = time.time()
        t_obj = delta_time(fprop) #get wanted properties in t=n

        '---- Get pressure field and new time step (if the past time step does \
        not obey the CFL condition) -------------------------------------------'

        self.delta_t = CompositionalFVM()(M, wells, fprop, self.delta_t, self.t)

        self.t += self.delta_t
        '----------------- Perform Phase stability test and flash -------------'

        if ctes.load_k and ctes.compressible_k:
            #if self.z

            self.p2 = StabilityCheck(fprop.P, fprop.T)
            fprop.L, fprop.V, fprop.xkj[0:ctes.Nc, 0, :], \
            fprop.xkj[0:ctes.Nc, 1, :], fprop.Csi_j[:,0,:], \
            fprop.Csi_j[:,1,:], fprop.rho_j[:,0,:], fprop.rho_j[:,1,:]  =  \
            self.p2.run_init(fprop.P, np.copy(fprop.z))#, wells)
            if len(wells['ws_inj'])>0 and not self.p2.constant_K:

                if any(wells['inj_cond']=='reservoir'):
                    injP_bool = wells['ws_p'] == wells['ws_inj']
                    injP = wells['ws_p'][injP_bool]

                    z = (wells['z']).T

                    p_well = StabilityCheck(fprop.P[wells['ws_inj']], fprop.T)
                    L, V, x, y, Csi_L, Csi_V, rho_L, rho_V  =  \
                    p_well.run_init(fprop.P[wells['ws_inj']],z[0:ctes.Nc])

                    if any(injP_bool):
                        qk_molar = fprop.qk_molar[:,injP]
                        #wells['values_q_vol'] = np.zeros_like((ctes.n_components,len(injP)))
                        wells['values_q_vol'] = qk_molar / \
                            ((Csi_V * V + Csi_L * L))#[wells['inj_cond']=='reservoir']
                    else:
                        wells['values_q'][:,wells['inj_cond']=='reservoir'] = (Csi_V * V + Csi_L * L) * self.q_vol
                    #wells['values_q_vol'][:,wells['inj_cond']=='reservoir'] = self.q_vol
                #else:
                    #wells['inj_p_term'] = []
                    #qk_molar = wells['values_q'][:,wells['inj_cond']=='surface']
                    #wells['values_q_vol'][:,wells['inj_cond']=='surface'] = qk_molar / \
                    #    ((Csi_V * V + Csi_L * L))[wells['inj_cond']=='surface']

        #if any(fprop.L!=0): import pdb; pdb.set_trace()
        '----------------------- Update fluid properties ----------------------'
        self.p1.run_inside_loop(M, fprop)
        #if (fprop.Sg[-37]<fprop.Sg[-36])*(fprop.Sg[-37]<fprop.Sg[-38]): import pdb; pdb.set_trace()


        '-------------------- Advance in time and save results ----------------'
        self.prod_rate_SC(fprop, wells)
        self.update_vpi(fprop, wells)

        self.update_loop()
        t1 = time.time()
        dt = t1 - t0
        self.sim_time += dt

        if self.use_vpi:
            if np.round(self.vpi,3) in self.vpi_save:
                self.update_current_compositional_results(M, wells, fprop) #ver quem vou salvar
        else:
            if self.time_save[0] == 0.0 or self.t in self.time_save:
                self.update_current_compositional_results(M, wells, fprop)

        self.delta_t = t_obj.update_delta_t(self.delta_t, fprop, wells, ctes.load_k, self.loop)#get delta_t with properties in t=n and t=n+1
        if len(wells['ws_prod'])>0: self.update_production(fprop, wells)


    def prod_rate_SC(self, fprop, wells):
        if ctes.load_k:
            p_well = StabilityCheck(ctes.P_SC, ctes.T_SC)
            z_prod = fprop.qk_prod/np.sum(fprop.qk_prod[0:ctes.Nc],axis=0)

            if (abs(z_prod.sum()-1)>1e-15) or np.isnan(z_prod.sum()): z_prod = fprop.z[:,wells['ws_prod']]
            L, V, x, y, Csi_L, Csi_V, rho_L, rho_V  =  \
                p_well.run_init(ctes.P_SC, z_prod[0:ctes.Nc])
            q_molar_prod = np.sum(fprop.qk_prod[0:ctes.Nc],axis=0)
            self.oil_production_rate_t = abs(np.sum(L * q_molar_prod / Csi_L))
            self.gas_production_rate_t = abs(np.sum(V * q_molar_prod / Csi_V))
        else: self.oil_production_rate_t =[]; self.gas_production_rate_t =[]

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
        return [np.array(['loop', 'vpi [s]', 'simulation_time [s]', 't [s]', \
                        'pressure [Pa]', 'Sw', 'So', 'Sg', 'Oil_p', 'Gas_p', \
                        'z', 'centroids', 'Nk', 'xkj', 'Oil production rate', \
                        'Gas production rate', 'time-step vector'])]

    def update_production(self, fprop, wells):
        ''' Function to compute oil and gas production [m³] through time'''

        if ctes.load_k:
            ws_p_prod = np.argwhere(wells['ws_p']==wells['ws_prod']).flatten()
            self.oil_production +=  abs(self.oil_production_rate_t)*self.delta_t
            self.gas_production +=  abs(self.gas_production_rate_t)*self.delta_t

            self.oil_production_rate.append(self.oil_production_rate_t)
            self.gas_production_rate.append(self.gas_production_rate_t)
            self.vector_time.append(self.t)

    def update_current_compositional_results(self, M, wells, fprop):
        #total_flux_internal_faces = fprop.total_flux_internal_faces.ravel() #* M.faces.normal[M.faces.internal]
        #total_flux_internal_faces_vector = fprop.total_flux_internal_faces.T * np.abs(M.faces.normal[M.faces.internal])
        if ctes.FR: Nk = fprop.Nk_SP;

        else: Nk = fprop.Nk

        self.current_compositional_results = np.array([self.loop, self.vpi, self.sim_time,
            self.t, fprop.P, fprop.Sw, fprop.So, fprop.Sg, self.oil_production,
            self.gas_production, fprop.z, M.data['centroid_volumes'], Nk, fprop.xkj,
            self.oil_production_rate, self.gas_production_rate, self.vector_time],dtype=object)
        self.all_results.append(self.current_compositional_results)
        M.data['saturation'][:] = fprop.Sw
        M.data['So'][:] = fprop.So
        M.data['pressure'][:] = fprop.P
        M.data['Sg'][:] = fprop.Sg
        #M.data['zC1'][:] = fprop.z[0,:]
        #M.data['perm'][:] = M.data[M.data.variables_impress['permeability']].reshape([ctes.n_volumes,3,3]).sum(axis=-1)[:,0]

        M.data.update_variables_to_mesh()
        # M.core.print(file = self.name_all_results + str(self.loop), extension ='.vtk')

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
