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
        
        if ctes.MUSCL['set']:
            from packs.compositional import prep_MUSCL as ctes_MUSCL
            ctes_MUSCL.run(M)

        fprop = self.get_initial_properties(M, wells)        

        return M, data_impress, wells, fprop, load, elements_lv0

    def run(self, M, wells, fprop, load, params, **kwargs):
        ''' Function created to compute reservoir and fluid properties at each \
        time step '''
        V_res = data_loaded['compositional_data']['water_data']['V_res']
        So_ini = data_loaded['compositional_data']['water_data']['So_ini']
        poro = data_loaded['compositional_data']['water_data']['poro']
        self.HCVP = V_res * fprop.So[-1] * poro

        t0 = time.time()
        t_obj = delta_time(fprop) #get wanted properties in t=n
        params['level0_negative_composition'][:] = False

        '---- Get pressure field and new time step (if the past time step does \
        not obey the CFL condition) -------------------------------------------'

        # kwargs['delta_t'] = self.delta_t
        params['delta_t'] = self.delta_t
        # self.delta_t = CompositionalFvmADM()(M, wells, fprop, **kwargs)
        self.delta_t = CompositionalFvmADM()(M, wells, fprop, delta_t=self.delta_t, t=self.t, params=params, **kwargs)

        self.t += self.delta_t
        
        print(f'\n Delta t: {self.delta_t} \n')

        '----------------- Perform Phase stability test and flash -------------'

        if ctes.load_k and ctes.compressible_k:
            #if self.z

            if ctes.load_k and ctes.compressible_k:
                self.p2 = StabilityCheck(fprop.P, fprop.T)

                fprop.L, fprop.V, fprop.A, fprop.xkj, fprop.Csi_j, fprop.rho_j = \
                    self.p2.run_init(fprop.P, np.copy(fprop.z), wells, \
                                     ksi_W=fprop.Csi_W, rho_W=fprop.rho_W)

                if ctes.miscible_w:
                    fprop.Csi_W = fprop.Csi_j[0, -1, :]
                    fprop.Csi_W0 = fprop.Csi_j[0, -1, :]
                    fprop.rho_W = fprop.Csi_j[0, -1, :]

            self.update_well_inj_rate(fprop, wells)


        '----------------------- Update fluid properties ----------------------'
        self.p1.run_inside_loop(M, fprop)


        '-------------------- Advance in time and save results ----------------'
        self.oil_production_rate_SC_t, self.gas_production_rate_SC_t, self.pressure_ = \
            self.prod_rate_RCorSC(fprop, wells, self.P_SC)
        self.oil_production_rate_RC_t, self.gas_production_rate_RC_t, self.pressure_ = \
            self.prod_rate_RCorSC(fprop, wells, fprop.P[wells['ws_prod']])

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
        if len(wells['ws_prod'])>0: self.update_production()
        '''Condições da WAG'''
        # import pdb; pdb.set_trace()
        if any(self.t >= wells['DT']):
            # import pdb; pdb.set_trace()
            if wells['z'][0][-1] == 0:
                wells['z'] = wells['za']
                wells['values_q'] = wells['values_qa']
                # self.get_well_inj_properties(M, fprop, wells)
            else:
                # import pdb; pdb.set_trace()
                wells['z'] = wells['zco2']
                wells['values_q'] = wells['values_qco2']
            # import pdb; pdb.set_trace()
            self.get_well_inj_properties(M, fprop, wells)
            wells['DT'] = wells['DT'][1:]

            # import pdb; pdb.set_trace()
