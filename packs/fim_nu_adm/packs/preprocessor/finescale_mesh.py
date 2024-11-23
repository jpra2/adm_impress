from .. import inputs
import numpy as np
# import line_profiler

def finescale_preprocess():
    finescale_inputs = inputs.finescale_inputs
    block_volume=np.array(finescale_inputs['mesh_generation_parameters']['block_size']).prod()
    GID_0, volume_centroids, face_areas, face_adjacencies = generate_centroids_areas_and_adjacencies(finescale_inputs['mesh_generation_parameters'])
    volume_permeabilities, porosities = get_permeability_and_porosity(volume_centroids, finescale_inputs['rock_parameters'], (finescale_inputs['mesh_generation_parameters']))
    face_permeabilities = get_k_harm_faces(volume_centroids, face_adjacencies, face_areas, volume_permeabilities)
    wells=get_wells(finescale_inputs['wells'], volume_centroids, GID_0)
    volumes={'GID_0':GID_0,'pore_volume':block_volume*porosities, 'Kxx':volume_permeabilities[:,0][GID_0],
             'centroids':volume_centroids}
    faces={'adjacent_volumes':face_adjacencies,'permeabilities':face_permeabilities}
    return wells, faces, volumes

def generate_centroids_areas_and_adjacencies(mesh):
    nb=mesh['n_blocks']
    lb=mesh['block_size']
    sp=mesh['starting_point']
    ms = np.mgrid[sp[0]+0.5*lb[0]:sp[0]+(nb[0]+0.5)*lb[0]:lb[0],
                  sp[1]+0.5*lb[1]:sp[1]+(nb[1]+0.5)*lb[1]:lb[1],
                  sp[2]+0.5*lb[2]:sp[2]+(nb[2]+0.5)*lb[2]:lb[2]]

    ms = ms.flatten()
    centroids = ms.reshape(3,int(ms.size/3)).T
    ############
    # ijk0 = np.array([centroids[:, 0]//lb[0], centroids[:, 1]//lb[1], centroids[:, 2]//lb[2]])
    # ee = ijk0[0] + ijk0[1]*nb[0] + ijk0[2]*nb[0]*nb[1] #60 3600
    # # import pdb; pdb.set_trace()
    # centroids=centroids[ee.astype(int)]
    # ##########
    GID_0 = np.arange(len(centroids))
    adjs0 = np.tile(GID_0,(np.array(nb)>1).sum())
    adjs1 = np.array([])
    sz=1
    sy=1
    if nb[2]>1:
        adjs1 = np.concatenate([adjs1, GID_0+1])
        sz=nb[2]
    if nb[1]>1:
        adjs1 = np.concatenate([adjs1, GID_0+sz])
        sy=nb[1]
    if nb[0]>1:
        adjs1 = np.concatenate([adjs1, GID_0+sz*sy])
    adjs1=adjs1.astype(int)
    adjs=np.vstack([adjs0,adjs1]).T
    adjs=adjs[adjs.max(axis=1)<=GID_0.max()]
    dif=abs(centroids[adjs[:,0]]-centroids[adjs[:,1]])
    adjacencies=adjs[((dif>0).sum(axis=1)==1) & ((dif<=lb).sum(axis=1)==3)]
    areas=np.array([lb[1]*lb[2], lb[0]*lb[2], lb[0]*lb[1]])
    areas=np.tile(areas,len(adjacencies)).reshape(len(adjacencies),3)[abs(centroids[adjacencies[:,0]]-centroids[adjacencies[:,1]])>0]
    return GID_0, centroids, areas, adjacencies

def get_permeability_and_porosity(centroids, rock, mgp):
    # import pdb; pdb.set_trace()
    if rock['constant_porosity'] or rock['constant_permeability']:
        block_dim=mgp['block_size']
        file_infos=rock['file_infos']
        data_spe10 = np.load(file_infos['name'])
        ks = data_spe10['perms'][:,[0,4,8]]
        phi = data_spe10['phi']
        phi = phi.flatten()
        nx , ny, nz = mgp['n_blocks']

        ijk0 = np.array([centroids[:, 0]//block_dim[0], centroids[:, 1]//block_dim[1], centroids[:, 2]//block_dim[2]])
        ee = ijk0[0] + ijk0[1]*60 + ijk0[2]*60*220 #60 3600
        # import pdb; pdb.set_trace()
        ee = ee.astype(np.int32)
        # import pdb; pdb.set_trace()
        permeabilities=ks[ee]
        porosities=phi[ee]

    if rock['constant_porosity']:
        porosities=np.repeat(rock['porosity_value'],len(centroids))
    if rock['constant_permeability']:
        permeabilities=np.tile(rock['permeability_value'],len(centroids)).reshape(centroids.shape)
        for area in rock["special_areas"]:
            region=rock["special_areas"][area]
            [rho,phi] = cart2pol((centroids[:,region['plane']]-np.array(region['center'])[region['plane']]).T)
            r0=region['radius']-region['thickness']/2
            r1=region['radius']+region['thickness']/2
            phi0, phi1=region['phi']
            pi=np.pi
            special_vols=(r1>rho)&(r0<rho)&(phi1*pi>phi)&(phi0*pi<phi)
            permeabilities[special_vols]=region['permeability']

    return permeabilities, porosities

def cart2pol(coords):
    x, y=coords
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def get_k_harm_faces(centroids, adjs, areas, ks):
    dist=centroids[adjs[:,1]]-centroids[adjs[:,0]]
    aux=dist>0
    l=dist.max(axis=1)
    k0 = ks[adjs[:,0]][aux]
    k1 = ks[adjs[:,1]][aux]
    k_harm = 2*k0*k1/(k0*l+k1*l)
    return k_harm

def get_wells(all_wells, centroids, GID_0):
    wells={}
    injectors=[]
    producers=[]
    diric=[]
    val_diric=[]
    neum=[]
    val_neum=[]
    bs=inputs.finescale_inputs['mesh_generation_parameters']['block_size']
    for w in all_wells:
        well=all_wells[w]
        p0=np.array(well['p0'])*bs
        p1=np.array(well['p1'])*bs

        vols = GID_0[((centroids>p0) & (centroids<p1)).sum(axis=1)==3]
        if len(vols>0):
            if well['type'] == 'Injector':
                injectors.append(vols)
            elif well['type'] == 'Producer':
                producers.append(vols)

            if well['prescription'] == 'P':
                diric.append(vols)
                val_diric.append(np.repeat(well['value'],len(vols)))
            elif well['prescription'] == 'Q':
                neum.append(vols)
                val_neum.append(well['value'])

    wells['ws_inj'] = np.concatenate(injectors)
    wells['ws_prod'] = np.concatenate(producers)
    wells['ws_p'] = np.concatenate(diric)
    wells['values_p'] = np.concatenate(val_diric)
    if len(neum)>0:
        wells['ws_q'] = np.concatenate(neum)
        wells['values_q'] = np.concatenate(val_neum)
    else:
        wells['ws_q'] = np.array([])
        wells['values_q'] = np.array([])
    return wells
