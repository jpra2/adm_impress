

def multilevel_pressure_solver(
        fine_scale_transmissibility,
        fine_scale_source_term,
        prolongation_list: list,
        restriction_list: list
):

    from packs.solvers.solvers_scipy.solver_sp import SolverSp
    solver = SolverSp()

    assert isinstance(prolongation_list, list)
    assert isinstance(restriction_list, list)
    n = len(prolongation_list)
    assert n == len(restriction_list)

    t2 = fine_scale_transmissibility.copy()
    q2 = fine_scale_source_term.copy()

    for OP, OR in zip(prolongation_list, restriction_list):
        t2 = t2 * OP
        t2 = OR * t2
        q2 = OR * q2

    solution = solver.direct_solver(t2, q2)
    prolongation_list.reverse()

    for OP in prolongation_list:
        solution = OP * solution

    return solution
