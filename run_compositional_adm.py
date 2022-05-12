from packs.directories import data_loaded
import pdb
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.IMPEC.compositionalIMPEC_adm import CompositionalFvmADM
from packs.compositional.update_time import delta_time
from packs.utils import constants as ctes
from run_compositional import run_simulation
import numpy as np
import time
# from packs.utils.test_functions import test_kwargs_keys, test_instance

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

class RunSimulationAdm(run_simulation):
    
    _kwargs_keys = {
        'run': set(['multilevel_data',
                     'multilevel_operators']),
    }

    def initialize(self, load, convert, mesh, **kwargs):
        ''' Function to initialize mesh (preprocess) get and compute initial mesh \
        properties '''
        M, elements_lv0, data_impress, wells = initial_mesh(mesh, load=load, convert=convert)
        ctes.init(M, wells)
        ctes.component_properties()
        if ctes.FR:
            from packs.compositional import prep_FR as ctes_FR
            ctes_FR.run(M)
        fprop = self.get_initial_properties(M, wells)        

        return M, data_impress, wells, fprop, load, elements_lv0

    def run(self, M, wells, fprop, load, params, **kwargs):
        ''' Function created to compute reservoir and fluid properties at each \
        time step '''

        t0 = time.time()
        t_obj = delta_time(fprop) #get wanted properties in t=n

        '---- Get pressure field and new time step (if the past time step does \
        not obey the CFL condition) -------------------------------------------'

        # kwargs['delta_t'] = self.delta_t
        params['delta_t'] = self.delta_t
        # self.delta_t = CompositionalFvmADM()(M, wells, fprop, **kwargs)
        self.delta_t = CompositionalFvmADM()(M, wells, fprop, delta_t=self.delta_t, t=self.t, params=params, **kwargs)

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