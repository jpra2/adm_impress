file='./mesh/4x4x4.msh'

# @profile
def importe():
    from impress import FineScaleMesh as mesh
    from packs.stokes_brinkman_3d import stokes_brinkman as stokes
    import multiprocessing as mp
    import gc
    import numpy as np
    import time
    # @profile
    def get_sol2(file):
        # gc.collect(generation=0)
        M = mesh(file)
        stokes.stokes_solver(M)
        # del(M)
    # @profile
    def run_process():
        p=mp.Process(target=get_sol2, args=[file])
        p.start()
        p.join()

    # run_process()
    get_sol2(file)
# import update_inputs
importe()
#
# @profile
# def get_sol(file):
#     M=mesh(file)
#     stokes.stokes_solver(M)
#     # gc.collect(generation=0)
#     # gc.is_tracked(M)
#
#
#     # gc.collect(generation=2)
# @profile
# def nada():
#     a=np.ones((1000,1000))
#     del(a)
#
# @profile
# def run_process():
#     p=mp.Process(target=get_sol2, args=[file])
#     p.start()
#     p.join()

    # import pdb; pdb.set_trace()




# run_process()
# nada()
# subDomains = [file]
# p=mp.Process(target=acumulate_OP, args=[file])
# p.start()
# p.join()
# procs = [mp.Process(target=acumulate_OP, args=[s]) for s in subDomains]


# get_sol(file)
# get_sol2(file)
# import pdb; pdb.set_trace()


# for i in range(5):
#     get_sol(file)



# import pdb; pdb.set_trace()
