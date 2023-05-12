from packs.manager.arraydatamanager import SuperArrayManager

bc_data = ['nodes_pressure_defined', 'nodes_pressure_values']


class BoundaryConditions(SuperArrayManager):
    __slots__ = bc_data