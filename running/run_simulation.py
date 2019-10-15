

class RunSimulation:

    def __init__(self, state=0, last=2):
        '''
        estado 0: gerar M
        estado 1: carregar M do estado 0
        estado 2: gerar malha dual
        estado 3: carregar malha dual
        '''

        self._loaded = False
        self.state = state
        steps = [0, 1, 2, 3]
        if state not in steps:
            raise ValueError(f'\nO estado: {state} nao esta em steps {steps}\n')

        while self.state <= last:
            if self._loaded:
                if self.state in [1, 3]:
                    # a malha esta no mesmo estado que foi gerada no estado 0
                    pass
                elif self.state in [2]:
                    self.run_dual_primal(M, self.state)

            else:
                if self.state in [0, 1]:
                    M = self.run_generate(self.state)
                elif self.state in [2, 3]:
                    M = self.run_generate(1)
                    self.run_dual_primal(M, self.state)

                    import pdb; pdb.set_trace()



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

    def step0_generate(self):
        from .run0 import M
        return M

    def step1_load(self):
        from .run1 import M
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

    def run_dual_primal(self, M, state):

        if state == 2:
            from .run2 import init_dual_mesh
            init_dual_mesh(M)
        elif state == 3:
            # TODO: carregar malha dual

            pass

    def run_load(self, state):
        pass
