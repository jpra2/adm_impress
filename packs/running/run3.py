from packs.multiscale import DualPrimalMesh1

def init_dual_mesh(M):
    dual_primal = DualPrimalMesh1(M, carregar=True)
    dual_primal.load_tags(M)
    dual_primal.get_elements(M)
