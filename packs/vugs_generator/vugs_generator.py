import random
class GenarateVugs:
    def __init__(self,M):
        centro=[300,500,700]
        abc=[1,3,2]
        self.get_comparation_parameters(M,centro,abc)
        centers, abcs=self.get_aleatory_parameters(10,[5,15]) #(num_vugs,[vmin, vmax])
        for i in range(len(centers)):
            self.get_comparation_parameters(M,centers[i],abcs[i])

        m=M.core.mb.create_meshset()
        M.core.mb.add_entities(m,M.core.all_volumes)
        M.core.mb.write_file("results/vugs.vtk",[m])


    def get_comparation_parameters(self,M, centro, params):
        centroids=M.volumes.center[:]
        x=centroids[:,0]
        y=centroids[:,1]
        z=centroids[:,2]
        self.x=x
        self.y=y
        self.z=z
        vugs=(x-centro[0])*(x-centro[0])/(params[0]*params[0])\
            +(y-centro[1])*(y-centro[1])/(params[1]*params[1])\
            +(z-centro[2])*(z-centro[2])/(params[2]*params[2])<1
        M.vug[vugs]=1
        
    def get_aleatory_parameters(self,n,lim_abc):
        xmin=self.x.min()
        xmax=self.x.max()
        ymin=self.y.min()
        ymax=self.y.max()
        zmin=self.z.min()
        zmax=self.z.max()
        centers=[]
        parameters=[]
        for i in range(n):
            centers.append([random.uniform(xmin, xmax),random.uniform(ymin,ymax), random.uniform(zmin,zmax)])
            parameters.append([random.uniform(lim_abc[0],lim_abc[1]),random.uniform(lim_abc[0],lim_abc[1]),random.uniform(lim_abc[0],lim_abc[1])])
        return centers,parameters
