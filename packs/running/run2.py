# criar malha dual

from packs.multiscale import DualPrimalMesh1

def init_dual_mesh(M):
    dual_primal = DualPrimalMesh1(M)
    dual_primal.run(M)

    pass
