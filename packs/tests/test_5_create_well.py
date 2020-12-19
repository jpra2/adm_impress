from packs.well_model.create_well import Well, AllWells
import numpy as np

volumes_ids = np.arange(10)
centroids = np.ones((len(volumes_ids) * 3)).reshape((len(volumes_ids), 3))
block_dimension = np.ones((len(volumes_ids), 3))
well = Well()
well.set_well(volumes_ids, centroids, block_dimension.copy(), direction='z')
well2 = Well()
volumes_ids2 = np.arange(23, 30)
block_dimension = np.ones((len(volumes_ids2), 3))
centroids = np.zeros((len(volumes_ids2), 3))
well2.set_well(volumes_ids2, centroids, block_dimension.copy())
permeability1 = np.tile(np.array([1, 0, 0, 0, 1, 0, 0, 0, 1]), len(volumes_ids)).reshape((len(volumes_ids), 3, -1))
well.update_req(well.calculate_req(well_permeability=permeability1))
well.wi = well.calculate_WI(well.blo)
well.update_n_phases(2)
well.phases = ['water', 'oil']
well.rho = [1000, 100]

import pdb;

pdb.set_trace()

AllWells.create_database()

# AllWells.load_database()
# all_wells = AllWells.load_wells_from_database()
