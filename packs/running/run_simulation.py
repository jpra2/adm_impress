from .. import directories as direc


class RunSimulation:

    def __init__(self, state=0, last=5):
        '''
        estado 0: gerar M
        estado 1: carregar M do estado 0
        estado 2: gerar malha dual
        estado 3: carregar malha dual
        estado 4: setar pocos
        estado 5: carregar pocos
        '''

        self._loaded = False
        self.state = state
        steps = [0, 1, 2, 3, 4, 5]
        if state not in steps:
            raise ValueError(f'\nO estado: {state} nao esta em steps {steps}\n')

        while self.state <= last:
            if self._loaded:
                if self.state in [1, 3, 5]:
                    # a malha esta no mesmo estado que foi gerada no estado atual
                    pass
                elif self.state in [2]:
                    self.step2_run_dual_primal(M, self.state)
                elif self.state in [4]:
                    self.step4_set_wells(M)

            else:
                if self.state in [0, 1]:
                    M = self.run_generate(self.state)
                elif self.state in [2]:
                    M = self.run_generate(1)
                    self.step2_run_dual_primal(M, self.state)
                elif self.state in [3]:
                    name_mesh = direc.names_outfiles_steps[2]
                    M = self.step1_load(name_mesh=name_mesh)
                    self.step3_load_dual(M)
                elif self.state in [4]:
                    name_mesh = direc.names_outfiles_steps[2]
                    M = self.step1_load(name_mesh=name_mesh)
                    self.step3_load_dual(M)
                    self.step4_set_wells(M)
                elif self.state in [5]:
                    name_mesh = direc.names_outfiles_steps[4]
                    M = self.step1_load(name_mesh=name_mesh)
                    self.step3_load_dual(M)
                    self.step5_load_wells(M)



                # ##################
                # # TODO: substituir:
                # if self.state == 0:
                #     M = self.run_generate(self.state)
                # elif self.state == 1:
                #     M = self.run_generate(self.state)
                #
                # por:
                # if self.state in [0, 1]:
                #     M = self.run_generate(self.state)
                #
                # pois eh a mesma coisa. so esta assim para fins didaticos
                # para compreender a sequencia de passos
                # ##################
                self._loaded = True
            self.state += 1

        self.M = M

    def step0_generate(self):
        print('\n--------------- Gerando Malha Inicial ---------------\n')
        from .run0 import M
        return M

    def step1_load(self, name_mesh=None):

        print('\n--------------- Carregando Malha Inicial ---------------\n')
        from .run1 import load_mesh
        if name_mesh:
            M = load_mesh(name_mesh = name_mesh)
        else:
            M = load_mesh()

        return M

    def run_generate(self, state):
        '''
        Gerar o objeto M
        '''
        if state == 0:
            ## setar malha inicial
            M = self.step0_generate()
            return M
        elif state == 1:
            M = self.step1_load()
            return M

    def step2_run_dual_primal(self, M, state=2):
        print('\n--------------- Criando Malha Dual ---------------\n')

        if state == 2:
            from .run2 import init_dual_mesh
            init_dual_mesh(M)

    def step3_load_dual(self, M, state=3):

        print('\n--------------- Carregando Malha Dual ---------------\n')

        from .run3 import init_dual_mesh
        init_dual_mesh(M)

    def step4_set_wells(self, M):
        print('\n--------------- Setando pocos ---------------\n')
        from .run4 import init_contours
        init_contours(M)

    def step5_load_wells(self, M):
        print('\n--------------- Carregando pocos ---------------\n')
        from .run5 import load_contours
        load_contours(M)

    def run_load(self, state):
        pass
