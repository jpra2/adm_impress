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

if data_loaded['water_miscible']:
    from packs.compositional.stability_check_3ph import StabilityCheck
    from packs.compositional.IMPEC.properties_calculation_3ph import PropertiesCalc
else:
    if data_loaded['compositional_data']['component_data']['constant_K']:
        from packs.compositional.Kflash import StabilityCheck
    else:
        from packs.compositional.stability_check import StabilityCheck

class run_simulation:
    '''Class created to compute simulation properties at each simulation time'''
    def __init__(self, name_current, name_all):
        #name_all = name_all + data_loaded['solver'] + hiperbolic_method + ctes.RS + ctes.
        self.name_current_results = os.path.join(direc.flying, name_current + '.npy')
        self.name_all_results = os.path.join(direc.flying, name_all)
        self.loop = 0
        self.vpi = 0.0
        self.t = 0.0
        self.oil_production = 0.
        self.gas_production = 0.
        self.oil_production_rate = []
        self.gas_production_rate = []
        self.cum_oil_prod = []

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

        if ctes.MUSCL['set']:
            from packs.compositional import prep_MUSCL as ctes_MUSCL
            ctes_MUSCL.run(M)

        fprop = self.get_initial_properties(M, wells)
        return M, data_impress, wells, fprop, load

    def get_well_inj_properties(self, M, fprop, wells):

        z = (wells['z']).T
        if (wells['ws_p']!=wells['ws_inj']): wells['inj_p_term'] = []
        #q_wells = wells['ws_inj']#[wells['inj_cond']=='reservoir']
        if ~ctes.miscible_w: z=z[0:ctes.Nc]
        if ctes.load_k and any(z.flatten()>0):
            if any(wells['inj_cond']=='reservoir'):
                P = fprop.P[wells['ws_inj']]
            else:
                P = np.array([ctes.P_SC])

            p_well = StabilityCheck(P, ctes.T_SC)
            L, V, A, xkj_ws, Csi_j_ws, rho_j_ws =  \
                p_well.run_init(P, z, ksi_W = fprop.Csi_W[wells['ws_inj']], \
                rho_W = fprop.rho_W[wells['ws_inj']])
            #import pdb; pdb.set_trace()
        else:
            L=np.zeros_like(wells['ws_inj']); V=np.copy(L)
            A = np.ones_like(wells['ws_inj'])
            Csi_j_ws = np.zeros_like(fprop.Csi_j[...,wells['ws_inj']])
            Csi_j_ws[-1,-1,:] = fprop.Csi_W[wells['ws_inj']]
            rho_j_ws = np.zeros_like(fprop.rho_j[...,wells['ws_inj']])
            rho_j_ws[-1,-1,:] = fprop.rho_W[wells['ws_inj']]
            xkj_ws = np.zeros_like(fprop.xkj[...,wells['ws_inj']])
            xkj_ws[-1,1,:] = 1

        'if the injector well has a prescribed pressure condition'

        ws_p_inj_ind = np.argwhere(wells['ws_p']==wells['ws_inj'])
        ws_p_inj = wells['ws_p'][ws_p_inj_ind].flatten()
        ws_inj_p = np.argwhere(wells['ws_inj']==wells['ws_p']).flatten()

        fprop.z[..., ws_p_inj] = z[:,ws_inj_p]

        mobility_ratio_ws = np.empty((1,ctes.n_phases,len(ws_p_inj)))

        mobility_ratio_ws[:,ctes.n_phases-1,:] = (A/(L+V+A))[ws_p_inj] #?? CHANGE THIS!!

        mobility_ratio_ws[:,0] = (L/(L+V+A))[ws_p_inj] #?? CHANGE THIS!!
        mobility_ratio_ws[:,1] = (V/(L+V+A))[ws_p_inj] #?? CHANGE THIS!!

        wells['inj_p_term'] = xkj_ws * mobility_ratio_ws * Csi_j_ws

        'if the injector well has a prescribed flux condition'
        if len(wells['ws_q'])>0:
            self.q_vol = np.copy(wells['values_q'])
            #rever esse Csi_L*L + Csi_V*V
            wells['values_q'] = (Csi_j_ws[:,1] * V + Csi_j_ws[:,0] * L +\
                Csi_j_ws[:,ctes.n_phases-1]*A) * self.q_vol
            wells['values_q_vol'] = self.q_vol
        #import pdb; pdb.set_trace()

    def get_initial_properties(self, M, wells):
        ''' get initial fluid - oil, gas and water data and calculate initial \
        properties'''

        fprop = FluidProperties(M, wells) # load reservoir properties data and initialize other data
        if ctes.load_w and not ctes.miscible_w:
            fprop.inputs_water_properties(M) #load water properties
        else:
            fprop.Csi_W = fprop.Csi_j[0,-1,:]
            fprop.Csi_W0 = fprop.Csi_j[0,-1,:]
            fprop.rho_W = fprop.Csi_j[0,-1,:]
        '------------------------- Perform initial flash ----------------------'
        self.get_well_inj_properties(M, fprop, wells)
        #fprop.z[:,0] = np.array([0,1])
        if ctes.load_k:
            self.p2 = StabilityCheck(fprop.P, fprop.T)
            fprop.L, fprop.V, fprop.A, fprop.xkj, fprop.Csi_j, fprop.rho_j =  \
            self.p2.run_init(fprop.P, np.copy(fprop.z), ksi_W = fprop.Csi_W, rho_W = fprop.rho_W)

            if ctes.miscible_w:
                fprop.Csi_W = fprop.Csi_j[0,-1,:]
                fprop.Csi_W0 = fprop.Csi_j[0,-1,:]
                fprop.rho_W = fprop.Csi_j[0,-1,:]

            ctes.P_SC *= np.ones_like(wells['ws_prod'])
            p_well = StabilityCheck(ctes.P_SC, ctes.T_SC)
            L, V, A, xkj, Csi_j, rho_j =  \
                p_well.run_init(ctes.P_SC, fprop.z[0:ctes.Nc,wells['ws_prod']], \
                ksi_W = fprop.Csi_W[wells['ws_prod']], rho_W = fprop.rho_W[wells['ws_prod']])

        else: fprop.xkj = []; fprop.L = []; fprop.V = []; wells['inj_term'] = []


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
        #if fprop.Sg[0]<1: import pdb; pdb.set_trace()

        if ctes.load_k and ctes.compressible_k:
            #if self.z
            self.p2 = StabilityCheck(fprop.P, fprop.T)

            fprop.L, fprop.V, fprop.A, fprop.xkj, fprop.Csi_j, fprop.rho_j =  \
            self.p2.run_init(fprop.P, np.copy(fprop.z), wells, \
            ksi_W = fprop.Csi_W, rho_W = fprop.rho_W)

            if ctes.miscible_w:
                fprop.Csi_W = fprop.Csi_j[0,-1,:]
                fprop.Csi_W0 = fprop.Csi_j[0,-1,:]
                fprop.rho_W = fprop.Csi_j[0,-1,:]

        self.update_well_inj_rate(fprop, wells)

        '----------------------- Update fluid properties ----------------------'
        self.p1.run_inside_loop(M, fprop)


        '-------------------- Advance in time and save results ----------------'
        self.prod_rate_SC(fprop, wells)
        self.update_vpi(fprop, wells)

        self.update_loop()
        t1 = time.time()
        dt = t1 - t0
        self.sim_time += dt
        #import pdb; pdb.set_trace()
        if self.use_vpi:
            if np.round(self.vpi,2) in self.vpi_save:
                self.update_current_compositional_results(M, wells, fprop) #ver quem vou salvar
        else:
            if self.time_save[0] == 0.0 or self.t in self.time_save:
                self.update_current_compositional_results(M, wells, fprop)

        self.delta_t = t_obj.update_delta_t(self.delta_t, fprop, wells, ctes.load_k, self.loop)#get delta_t with properties in t=n and t=n+1
        if len(wells['ws_prod'])>0: self.update_production(fprop, wells)
        #import pdb; pdb.set_trace()

    def update_well_inj_rate(self, fprop, wells):
        # Update wells injection rate if injection condition is the reservoir one
        ws_inj_reservoir = wells['inj_cond']=='reservoir'
        if any(ws_inj_reservoir) and not self.p2.constant_K:
            # look for injection wells with prescribed pressure
            injP_bool = wells['ws_p'] == wells['ws_inj']
            injP = wells['ws_p'][injP_bool]

            z = (wells['z']).T
            if not ctes.miscible_w: z = z[0:ctes.Nc]
            if any(z>0):
                p_well = StabilityCheck(fprop.P[wells['ws_inj']], fprop.T)
                L, V, A, xkj, Csi_j, rho_j  =  \
                p_well.run_init(fprop.P[wells['ws_inj']],z, \
                    ksi_W = fprop.Csi_W[wells['ws_inj']],
                    rho_W = fprop.rho_W[wells['ws_inj']])

                if any(injP_bool): #tentar remover esse if!
                    qk_molar = fprop.qk_molar[:,injP]
                    #wells['values_q_vol'] = np.zeros_like((ctes.n_components,len(injP)))
                    wells['values_q_vol'] = qk_molar / \
                        ((Csi_j[:,1] * V + Csi_j[:,0] * L + Csi_j[:,ctes.n_phases-1] * A)[:,np.newaxis])#[wells['inj_cond']=='reservoir']
                else:
                    wells['values_q'][:,wells['inj_cond']=='reservoir'] = (Csi_j[:,1] * V + \
                        Csi_j[:,0] * L+ Csi_j[:,ctes.n_phases-1]*A) * self.q_vol
                #wells['values_q_vol'][:,wells['inj_cond']=='reservoir'] = self.q_vol
            else:
                wells['values_q'][:,wells['inj_cond']=='reservoir'] = (fprop.Csi_j[:,-1,wells['inj_cond']=='reservoir']) * self.q_vol

    def prod_rate_SC(self, fprop, wells):
        if ctes.load_k:
            p_well = StabilityCheck(ctes.P_SC, ctes.T_SC)
            z_prod = fprop.qk_prod/np.sum(fprop.qk_prod[0:ctes.Nc],axis=0)
            if ~ctes.miscible_w: z_prod=z_prod[0:ctes.Nc]
            if (abs(z_prod.sum()-1)>1e-15) or np.isnan(z_prod.sum()): z_prod = fprop.z[:,wells['ws_prod']]
            L, V, A, xkj, Csi_j, rho_j  =  \
                p_well.run_init(ctes.P_SC, z_prod[0:ctes.Nc], \
                ksi_W = fprop.Csi_W[wells['ws_prod']], rho_W = fprop.rho_W[wells['ws_prod']])
            q_molar_prod = np.sum(fprop.qk_prod[0:ctes.Nc],axis=0)
            #import pdb; pdb.set_trace()
            self.oil_production_rate_t = abs(np.sum(L * q_molar_prod / Csi_j[:,0]))
            self.gas_production_rate_t = abs(np.sum(V * q_molar_prod / Csi_j[:,1]))
            if abs(q_molar_prod)>1e3: import pdb; pdb.set_trace()

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
        #flux_total_inj = abs(fprop.qt_inj)
        #import pdb; pdb.set_trace()
        #print(self.vpi)

        self.vpi = self.vpi + (flux_total_inj.sum())/sum(fprop.Vp)*self.delta_t
        #print(self.vpi)
        #if self.vpi>0.5: import pdb; pdb.set_trace()

    def get_empty_current_compositional_results(self):
        return [np.array(['loop', 'vpi [s]', 'simulation_time [s]', 't [s]', \
                        'pressure [Pa]', 'Sw', 'So', 'Sg', 'Oil_p', 'Gas_p', \
                        'z', 'centroids', 'Nk', 'xkj', 'Oil production rate', \
                        'Gas production rate', 'Cumulative oil production',
                         'time-step vector'])]

    def update_production(self, fprop, wells):
        ''' Function to compute oil and gas production [m³] through time'''

        if ctes.load_k:
            ws_p_prod = np.argwhere(wells['ws_p']==wells['ws_prod']).flatten()
            self.oil_production +=  abs(self.oil_production_rate_t)*self.delta_t
            self.gas_production +=  abs(self.gas_production_rate_t)*self.delta_t

            self.oil_production_rate.append(self.oil_production_rate_t)
            self.gas_production_rate.append(self.gas_production_rate_t)
            self.vector_time.append(self.t)
            self.cum_oil_prod.append(self.oil_production)

    def update_current_compositional_results(self, M, wells, fprop):
        #total_flux_internal_faces = fprop.total_flux_internal_faces.ravel() #* M.faces.normal[M.faces.internal]
        #total_flux_internal_faces_vector = fprop.total_flux_internal_faces.T * np.abs(M.faces.normal[M.faces.internal])
        if ctes.FR: Nk = fprop.Nk_SP;

        else: Nk = fprop.Nk

        self.current_compositional_results = np.array([self.loop, self.vpi, self.sim_time,
            self.t, fprop.P, fprop.Sw, fprop.So, fprop.Sg, self.oil_production,
            self.gas_production, fprop.z, M.data['centroid_volumes'], Nk, fprop.xkj,
            self.oil_production_rate, self.gas_production_rate, self.cum_oil_prod,
            self.vector_time],dtype=object)
        self.all_results.append(self.current_compositional_results)
        M.data['saturation'][:] = fprop.Sw
        M.data['So'][:] = fprop.So
        M.data['pressure'][:] = fprop.P
        M.data['Sg'][:] = fprop.Sg


        #M.data['zC1'][:] = fprop.z[0,:]
        #M.data['perm'][:] = M.data[M.data.variables_impress['permeability']].reshape([ctes.n_volumes,3,3]).sum(axis=-1)[:,0]

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
