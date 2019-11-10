
def load_contours(M):
    from ..contours.set_contours import Contours
    contours = Contours(M)
    contours.load_tags(M)
    contours.load_from_npz()
    contours.loaded()
