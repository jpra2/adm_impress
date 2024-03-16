from packs.manager.arraydatamanager import SuperArrayManager
from copy import deepcopy
from packs.errors import err
from packs import defnames
import numpy as np


class BoundaryConditions(SuperArrayManager):
    """Boundary conditions class

        Set the boundary conditions of problem.
        The data name must be in cls.boundary_names.
        The data is an structured array where it is composed of an 'id' key and 'value'  key
        
        example:
             bc = BoundaryContitions()
             bc.insert_name('problem_name')
             dtype_bc_array = [('id', np.uint64), ('value', np.float64)]
             bc_array = np.zeros(2, dtype=dtype_bc_array)
             bc_array['id'][:] = [0, 1] ## nodes ids
             bc_array['value'][:] = [100, 0] ## pressures values ('nodes_pressures' data)
             bc.insert_data({'nodes_pressures': bc_array})

    Raises:
        err.NameExistsError: _description_
    """
    
    boundary_names = defnames.mpfa_boundary_names
    dtype_bc_array = [('id', np.uint64), ('value', np.float64)]
    
    @classmethod
    def test_names(cls, names):
        
        names2 = deepcopy(names)
        try:
            names2.remove('name')
        except ValueError:
            pass
        
        for name in names2:
            if name not in cls.boundary_names:
                raise err.NameExistsError(f'the name {name} not in {cls.boundary_names}')

    def set_boundary(self, bc_type, ids, values):
        self.test_names([bc_type])
        
        bc_array = np.zeros(len(ids), dtype=self.dtype_bc_array)
        bc_array['id'][:] = ids
        bc_array['value'][:] = values
        
        self.insert_data({bc_type: bc_array})

    def get_neumann_nodes(self, nodes_of_edges):
        neumann_edges = self['neumann_edges']['id']
        neumann_nodes = np.unique(nodes_of_edges[neumann_edges])
        neumann_nodes = np.setdiff1d(neumann_nodes, self['dirichlet_nodes']['id'])
        return neumann_nodes

