# from packs.tests import test_dataset
# from packs.tests import test_geometric_data
# from packs.tests import test_preprocess_tpfa
# from packs.tests import test_biphasic_tpfa_cons
# from packs.tests import test_2_generate_non_uniform_mesh
# from packs.tests import test_3_prolongation_with_pcorr
# from packs.tests import test_4_monophasic
# from packs.examples import mpfad_gls_test
# from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import test_weights
# test_weights()

from packs.examples.remove_data import run_delete_all

# from packs.mpfa_methods.flux_calculation.test.test_flux_lsds_method import test_lsds_flux
# test_lsds_flux()

# from packs.mpfa_methods.test.test_monophasic_lsds_method import test_monophasic_problem_with_pressure_prescription
# test_monophasic_problem_with_pressure_prescription()

# from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import test_neumann_weights
# test_neumann_weights()

# from packs.mpfa_methods.test.test_monophasic_lsds_dong_paper import plot_errors
# plot_errors()

# from packs.mpfa_methods.test.test_inverse_distance_weights import test_1
# test_1()

# from packs.mpfa_methods.weight_interpolation.test.test_lpew_weights import test_weights_prof_fernando, sequence
# test_weights_prof_fernando()

# from packs.mpfa_methods.test.test_monophasic_lsds_dong_paper import plot_errors
# # list_problems = [0, 1, 2, 3, 4]
# list_problems = [0]
# for problem in list_problems:
#     print(f'Running problem: {problem} \n')
#     plot_errors(problem)
#     print(f'Finish problem: {problem} \n')

# from packs.mpfa_methods.flux_calculation.test.test_flux_diamond_method import test_xi_params_ds
# test_xi_params_ds()

# from packs.mpfa_methods.test.test1_lpew2_monophasic import test_problem1, test_problem2, test_problem3, test_problem4
# test_problem2()
# test_problem1()
# test_problem3()
# test_problem4()

# from packs.mpfa_methods.test.test2_lpew2_monophasic import run
# run()

# from packs.manager.test.meshiowrapper_test import run
# run()

# from packs.multiscale.unstructured.test.create_primal_test import run
# run()

# from packs.preprocess.test.mesh_properties_from_meshio_test import run
# msh = run()

# from packs.multiscale.unstructured.test.create_primal_test import run
# run()

# from packs.multiscale.unstructured.test.test_create_dual import run2
# run2()

# from packs.multiscale.unstructured.test.test_uns_ams_prolongation import run3
# ams, fine, coarse = run3()

# from packs.mpfa_methods.weight_interpolation.test.test_gls_edges_weights import run
# run()

# # run_delete_all()
# from packs.multiscale.unstructured.test.test_brazil import run4
# run4()

# run_delete_all()
# from packs.multiscale.unstructured.test.test_cross import run4
# run4()

# from packs.tests.test_initial_preprocess import run
# run()

# from packs.examples.remove_data import run_flying, run_results
# run_flying()
# run_results()
# import pdb; pdb.set_trace()

# from packs.examples.biphasic_mpfa import run
# run()

# from packs.examples.biphasic_mpfa_nu_adm import run
# run()
# plot_graph()

# from packs.examples.monophasic_unstructured_linear import run4
# run4()

# from packs.examples.diss_test1 import run4
# run4()

# from packs.examples.biphasic_mpfa_nuadm_v2 import run
# run()

# run_delete_all()
from packs.examples.biphasic_mpfa_layers_tri_finescale import run5
run5()


