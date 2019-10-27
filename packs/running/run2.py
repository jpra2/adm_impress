# criar malha dual

from ..multiscale.dual_primal.create_dual_and_primal_mesh import DualPrimalMesh1

def init_dual_mesh(M):
    dual_primal = DualPrimalMesh1(M)
    dual_primal.run(M)

    pass
