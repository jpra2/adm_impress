
class EmptyQueueError(Exception):
    pass

class DualStructureError(Exception):
    pass

class ConservativeVolumeError(Exception):
    pass

class PmsFluxFacesError(Exception):
    pass

class MaxLoopIterationError(Exception):
    pass

class NameLoadError(Exception):
    pass



class NameExistsError(Exception):
    pass

class TagNameExistsError(Exception):
    pass

class ElementTypeNotInMeshError(Exception):
    pass

class DimensionError(Exception):
    pass

class MoabTypeNotFoundError(Exception):
    pass

class MinMaxValueTestError(Exception):
    pass
