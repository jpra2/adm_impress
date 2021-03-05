
def init_contours(M):
    from ..contours.set_contours import Contours
    contours = Contours(M)
    contours.create_tags(M)
    contours.get_wells(M)
    contours.correct_wells()
    contours.set_infos(M)
    contours.export_to_npz()
    contours.save_mesh(M)
    contours.loaded()
