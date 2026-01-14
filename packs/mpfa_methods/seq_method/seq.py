import numpy as np
from typing import Any
import torch
from torch.func import vjp, vmap, jacrev, jvp
import networkx as nx
import pandas as pd
import time
from scipy.sparse import csr_matrix
from scipy.sparse import linalg

from packs.solvers.solvers_scipy.solver_sp import SolverSp

def residual_function(v1: torch.TensorType) -> torch.TensorType:
    f1 = torch.zeros_like(v1)
    n = v1.size()[0]
    f1[0] = v1[0] - 100
    f1[n-1] = v1[n-1] - 1
    for i in range(1, n-1):
        # f1[i] = -v1[i-1] + 2*v1[i] - v1[i+1]
        f1[i] = v1[i]**2 - v1[i+1] - v1[i-1]
         
    return f1

def get_adj_matrix(v1: torch.TensorType) -> np.ndarray:
    
    n = v1.size()[0]
    adj = np.zeros((n,n), dtype=np.int64)
    adj[0, [0, 1]] = 1
    adj[n-1, [n-1, n-2]] = 1
    rows = [0, 0, n-1, n-1]
    cols = [0, 1, n-2, n-1]
    for i in range(1, n-1):
        # adj[i, i+1] = 1
        # adj[i, i-1] = 1
        # adj[i, i] = 1
        adj[i, [i-1, i, i+1]] = 1
        rows.extend([i, i, i])
        cols.extend([i-1, i, i+1])
    
    return adj, rows, cols
    

def get_mask(v1: torch.TensorType) -> torch.TensorType:
    n = v1.size()[0]
    adj = torch.zeros((n, n), dtype=torch.bool)
    adj[0, 0] = True
    adj[n-1, n-1] = True
    for i in range(1, n-1):
        adj[i, [i-1, i, i+1]] = True
    
    return adj

def mount_jacobian_v2(R, x: torch.Tensor, mask: torch.Tensor):
    n = x.size()[0]
    def vjp_row(v):
        y = R(x)
        return torch.autograd.grad(
            outputs=y,
            inputs=x,
            grad_outputs=v,
            retain_graph=True
        )[0]
    
    E = torch.eye(n)
    J = vmap(vjp_row)(E)
    import pdb; pdb.set_trace()
    rows, cols = mask.nonzero(as_tuple=True)
    
    values = J[rows, cols]

    J_sparse = torch.sparse_coo_tensor(
        torch.stack([rows, cols]),
        values,
        size=(n, n)
    ).to_sparse_csr()
    
    return J_sparse
    

def mount_jacobian_v1(R: torch.Tensor, x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    n = x.size()[0]
    rows, cols = mask.nonzero(as_tuple=True)
    values = torch.zeros(len(rows), dtype=torch.float64)
    unique_rows = torch.unique(rows)
    y = R(x)
    local_mask = rows == 0
    for i in unique_rows:
        local_mask[:] = rows == i
        local_cols = cols[local_mask]
        grad = torch.autograd.grad(
            outputs=y[i],
            inputs=x,
            retain_graph=True,
            allow_unused=True            
        )[0]
        values[local_mask] = grad[local_cols]
    
    J_sparse = torch.sparse_coo_tensor(
        torch.stack([rows, cols]),
        values,
        size=(n, n)
    ).to_sparse_csr()
    
    return J_sparse

def solve_jacobian(J: torch.Tensor, x0: torch.Tensor, b: torch.Tensor) -> np.ndarray:
    crow_indices = J.crow_indices().numpy() # indptr in SciPy terminology
    col_indices = J.col_indices().numpy() # indices in SciPy terminology
    values = J.values().numpy()
    shape = J.size()
    
    A2 = csr_matrix((values, col_indices, crow_indices), shape=shape)
    
    x1 = x0.detach().numpy()
    b1 = b.detach().numpy()
    
    dx, exitcode = linalg.bicgstab(A=A2, b=b1, x0=x1, atol=1e-10)
    return dx
    

def loop_Newton_method(R, x: torch.Tensor, b, mask, maxit=1000, tol=1e-8) -> torch.Tensor:
    
    y = R(x)
    nm = torch.norm(y)
    loop = 0
    x0 = torch.zeros_like(x)   
    w = 2/3
    
    while nm > tol and loop < maxit:
        # J = mount_jacobian_v1(R, x, mask)
        J = mount_jacobian_v2(R, x, mask)
        dx = solve_jacobian(J, x0, -y)
        with torch.no_grad():
            x[:] +=  dx
        y[:] = R(x)
        nm = torch.norm(y)
        print(nm.item())
        loop += 1        
    
    return x
        
        
    
    
    

def run():
    n = 1000
    x = torch.arange(n, requires_grad=True, dtype=torch.float64)
    b = torch.zeros(n, dtype=torch.float64)
    b[0] = 100
    b[n-1] = 1
    
    adj_matrix = get_mask(x)
    
    x = loop_Newton_method(residual_function, x, b, adj_matrix)
    
    
    
    
    
    
    
    
    
    
    
    # # 3. Compute compressed Jacobian using Vector-Jacobian Products
    # _, vjp_fn = vjp(residual_function, v1)
    # compressed_jacobian = vmap(vjp_fn)(seed_matrix.T)
    # jac2 = jacrev(residual_function)(v1)
    
    # adj_matrix, rows, cols = get_adj_matrix(v1)
    # t0 = time.perf_counter()
    # adj_matrix = get_mask(v1)
    # rows, cols = adj_matrix.nonzero(as_tuple=True)
    # # G = nx.from_numpy_array(adj_matrix.numpy())
    # # colors = nx.coloring.greedy_color(G, strategy="largest_first")
    # # indices = torch.tensor([rows, cols], dtype=torch.int32)
    
    
    # # nodes = np.array(list(colors.keys()))
    # # nodes_colors = np.array(list(colors.values()))
    
    # # colormap = pd.Series(data=nodes_colors, index=nodes)
    # # colormap.sort_index(inplace=True)
    
    # # nodes = torch.from_numpy(colormap.index.values)
    # # nodes_colors = torch.from_numpy(colormap.values)
    
   
    
    # values = torch.zeros(len(rows))
    # unique_rows = torch.unique(rows)
    # y = residual_function(v1)
    
    # # for k, (i, j) in enumerate(zip(rows, cols)):
    # #     grad = torch.autograd.grad(
    # #         outputs=y[i],
    # #         inputs=v1,
    # #         retain_graph=True,
    # #         allow_unused=True            
    # #     )[0]
    # #     values[k] = grad[j] if grad is not None else 0.0
    
    
    # local_mask = rows == 0
    # for i in unique_rows:
    #     local_mask[:] = rows == i
    #     local_cols = cols[local_mask]
    #     grad = torch.autograd.grad(
    #         outputs=y[i],
    #         inputs=v1,
    #         retain_graph=True,
    #         allow_unused=True            
    #     )[0]
    #     values[local_mask] = grad[local_cols]
    
    # # for k, (i, j) in enumerate(zip(rows, cols)):
    # #     grad = torch.autograd.grad(
    # #         outputs=y[i],
    # #         inputs=v1,
    # #         retain_graph=True,
    # #         allow_unused=True            
    # #     )[0]
    # #     values[k] = grad[j] if grad is not None else 0.0
    
    # # rows2 = []
    # # cols2 = []
    # # vals2 = []
    
    # # for c in range(n_colors):
    # #     m1 = nodes_colors == c
    # #     v = m1.float()
    # #     _, Jv = jvp(residual_function, (v1,), (v,))
    # #     nz = Jv.nonzero(as_tuple=True)
        
    # #     # j = torch.where(v != 0)[0][0]
        
    # #     # for i in nz[0]:
    # #     #     j = torch.where(v != 0)[0][0]
    # #     #     rows2.append(i.item())
    # #     #     cols2.append(j.item())
    # #     #     vals2.append(Jv[i].item())
        
    # #     for k, (i, j) in enumerate(zip(rows, cols)):
    # #         if v[j] != 0:
    # #             values[k] = Jv[i]
    
    # J_sparse = torch.sparse_coo_tensor(
    #     torch.stack([rows, cols]),
    #     values,
    #     size=(n, n)
    # )
    # t1 = time.perf_counter()
    # el1 = t1 - t0
    # print(el1)
    
    
    
    
    
    
    
    # num_colors = torch.unique(nodes_colors).size()[0]
    # n = v1.size()[0]
    
    # # groups = []
    # # for c in torch.unique(nodes_colors):
    # #     test = nodes_colors == c
    # #     groups.append(nodes[test])
    
    # seed_matrix = torch.zeros((num_colors, n), dtype=torch.int64)
    # for c in range(num_colors):
    #     seed_matrix[c, nodes_colors == c] = 1.0
    
    # _, vjp_fn = vjp(residual_function, v1)
    # compressed_jac = vmap(vjp_fn)(seed_matrix)[0]
    
    # entry_colors = nodes_colors[indices[1]]
    # # v1 = entry_colors.int().detach().cpu().numpy()
    # # v2 = indices[0].int().detach().cpu().numpy()
    # sparse_vals = compressed_jac[entry_colors, indices[0]]    
    # sparse_jac = torch.sparse_coo_tensor(indices, sparse_vals, size=(n, n)).coalesce()
    
    # print(sparse_jac.to_dense())
    # print(jac2)
    # print(J_sparse.to_dense())
    
    
    
    
    
    # print(compressed_jacobian) 
    # print(jac2)
    # print(colors)
    
    
    
    
    
    
    