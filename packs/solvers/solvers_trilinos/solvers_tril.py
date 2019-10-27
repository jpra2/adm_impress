from PyTrilinos import Epetra, AztecOO, EpetraExt, Teuchos, IFPACK
import numpy as np
import scipy.sparse as sp
# IFPACK.PrintSparsity(Matrix, "matrix.ps")


class PyTrilWrap:

    def __init__(self, p=1):
        self._comm  = Epetra.PyComm()
        self._params = dict()
        self.set_parameters()

    def solve_linear_problem(self, A, b, x=None, its=1000, tolerance=1e-10):
        '''
        resolve o problema Ax = b
        input:
            A: matriz quadrada do scipy
            b = termo fonte
            x: chute inicial
            its: numero maximo de iteracoes
            tolerance: tolerancia para o residuo
        output:
            res: informa se o residuo foi menor que a tolerancia
            x: vetor resposta
        '''
        comm = self.comm
        n = len(b)
        std_map = Epetra.Map(n, 0, comm)
        x2 = Epetra.Vector(std_map)
        if x:
            x2[:] = x[:]
        b2 = Epetra.Vector(std_map)
        b2[:] = b[:]
        A2 = Epetra.CrsMatrix(Epetra.Copy, std_map, 7)
        indices = sp.find(A)
        A2.InsertGlobalValues(indices[0], indices[1], indices[2])
        irr = A2.FillComplete()
        linearProblem = Epetra.LinearProblem(A2, x2, b2)
        solver = AztecOO.AztecOO(linearProblem)
        solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_warnings)
        solver.SetParameters(self._params)
        solver.Iterate(its, tolerance)
        x2 = np.array(x2)
        res = solver.ScaledResidual() < tolerance
        return x2, res

    def set_parameters(self, params=None):
        if params:
            pass
        else:
            params = {'Solver': 'GMRES',
                      'Precond': 'Jacobi'}

        self._params.update(params)

'''
>>> solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_dom_decomp)
>>> solver.SetAztecOption(AztecOO.AZ_subdomain_solve, AztecOO.AZ_ilu)
>>> solver.SetAztecOption(AztecOO.AZ_overalp, 1)
>>> solver.SetAztecOption(AztecOO.AZ_graph_fill, 1)

>>> solver.SetParameters({"precond": "dom_decomp",
...                       "subdomain_solve": "ilu",
...                       "overlap": 1,
...                       "graph_fill": 1})

'''
