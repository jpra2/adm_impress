from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh

def load_m_object(mesh_name, dimension=3):
    M = msh(mesh_name, dim = dimension)
    print(f'\n{mesh_name} object loaded\n')
    return M

if __name__ == "__main__":
    mesh_name = '9x9x9.h5m'
    M = load_m_object(mesh_name, 3)
