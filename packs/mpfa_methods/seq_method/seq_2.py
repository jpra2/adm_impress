import numpy as np
from typing import Any
import torch
from torch.func import vjp, vmap, jacrev
import networkx as nx
import time

from .seq import residual_function, get_adj_matrix, get_mask


def run2():
    v1 = torch.arange(10000, requires_grad=True, dtype=torch.float32)
    v2 = torch.arange(5, requires_grad=True, dtype=torch.float32)
    n = v1.size()[0]
    
    with torch.no_grad():
        v2[2] = 5
        v2[3] = 10
    
    # f1 = torch.zeros_like(v1)
    # f1[0] = v1[0]**2 + v1[1]**2
    # f1[1] = v1[1]**2 + v1[2]**2
    # f1[2] = v1[2]**2 + v1[0]**2
    
    # groups = [torch.tensor([0, 1]), torch.tensor([1, 2]), torch.tensor([0, 2])]
    # seed_matrix = torch.zeros((3, len(groups)))
    # for i, group in enumerate(groups):
    #     seed_matrix[group, i] = 1.0
    
    # # 3. Compute compressed Jacobian using Vector-Jacobian Products
    # _, vjp_fn = vjp(residual_function, v1)
    # compressed_jacobian = vmap(vjp_fn)(seed_matrix.T)
    t0 = time.perf_counter()
    # jac2 = jacrev(residual_function)(v1)
    t1 = time.perf_counter()
    el1 = t1 - t0
    
    # adj_matrix, rows, cols = get_adj_matrix(v1)
    t0 = time.perf_counter()
    mask = get_mask(v1)
    
    def vjp_row(v):
        y = residual_function(v1)
        return torch.autograd.grad(
            outputs=y,
            inputs=v1,
            grad_outputs=v,
            retain_graph=True
        )[0]


    E = torch.eye(n)
    J = vmap(vjp_row)(E)
    
    rows, cols = mask.nonzero(as_tuple=True)
    
    values = J[rows, cols]

    J_sparse = torch.sparse_coo_tensor(
        torch.stack([rows, cols]),
        values,
        size=(n, n)
    )
    t1 = time.perf_counter()
    el2 = t1 - t0
    
    # print(f'Time traditional Jacobian: {el1}')
    print(f'Time masked Jacobian: {el2}')
    
    
    
    
    
    # print(compressed_jacobian) 
    # print(jac2)
    # print(colors)
    
    
    
    
    
    
    