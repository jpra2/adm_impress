from packs.preprocess.create_mesh import createMesh

mesh = createMesh()
mesh.create_fine_vertices()
mesh.create_elements()
mesh.export_mesh()
