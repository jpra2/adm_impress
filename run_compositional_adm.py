from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.IMPEC.compositionalIMPEC_adm import CompositionalFvmADM
from packs.compositional.update_time import delta_time
from packs.utils import constants as ctes
from run_compositional import run_simulation
import numpy as np
import time
# from packs.utils.test_functions import test_kwargs_keys, test_instance

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
        fprop = self.get_initial_properties(M, wells)

        return M, data_impress, wells, fprop, load, elements_lv0

    def run(self, M, wells, fprop, load, **kwargs):
        ''' Function created to compute reservoir and fluid properties at each \
        time step '''
        # adm_method = kwargs.get('adm_method')
        # neumann_subds = kwargs.get('neumann_subds')
        # data_impress = kwargs.get('data_impress')

        # test_kwargs_keys(self._kwargs_keys['run'], kwargs.keys())
        # params = kwargs.get('params')

        t0 = time.time()
        t_obj = delta_time(fprop) #get wanted properties in t=n

        '---- Get pressure field and new time step (if the past time step does \
        not obey the CFL condition) -------------------------------------------'

        # self.delta_t = CompositionalFvmADM()(M, wells, fprop, self.delta_t,
        #                                   multilevel_data=kwargs.get('multilevel_data'),
        #                                   multilevel_operators=kwargs.get('multilevel_operators'),
        #                                   params=params,
        #                                   adm_method=adm_method,
        #                                   neumann_subds=neumann_subds,
        #                                   data_impress=data_impress)

        # self.delta_t = CompositionalFvmADM()(M, wells, fprop, self.delta_t, **kwargs)
        kwargs['delta_t'] = self.delta_t
        self.delta_t = CompositionalFvmADM()(M, wells, fprop, **kwargs)

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