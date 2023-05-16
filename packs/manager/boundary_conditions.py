from packs.manager.arraydatamanager import SuperArrayManager

bc_data = ['nodes_pressures']


class BoundaryConditions(SuperArrayManager):
    __slots__ = bc_data