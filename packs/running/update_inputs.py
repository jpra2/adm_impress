from packs.utils.info_manager import InfoManager
import os
import pdb

def Update_variable(variable_file):

    # dd = InfoManager('input_cards/inputs0_2.yml', 'input_cards/inputs0.yml')
    # dd = InfoManager('input_cards/inputs0.yml', 'input_cards/inputs0.yml')
    dd = InfoManager(variable_file, 'input_cards/variable_input.yml')
    dd.save_obj()
    # dd['load_data'] = True
    # dd['multilevel_data'] = False
    # dd['load_multilevel_data'] = True
    # dd['load_operators'] = True
    # dd['load_biphasic_data'] = True
    # dd['capillary_pressure'] = True


    # dd['read_permeability'] = True
    # dd['set_permeability'] = False
    # dd['n_levels'] = 1

    # dd['biphasic'] = True

    # dd['deletar_results'] = True
    # dd['gravity'] = True
    # dd['get_correction_term'] = True
    # dd['convert_english_to_SI'] = True
    # dd.save_obj()

    # dd2 = InfoManager('input_cards/variable_adm.yml', 'input_cards/variable_input.yml')
    # dd2.save_obj()



def Delete_Files():
    dd = InfoManager('input_cards/inputs0.yml', 'input_cards/inputs0.yml')
    if dd['deletar_results'] == True:

        results = 'results'
        ff = os.listdir(results)

        for f in ff:
            # if f[-4:] == '.vtk':
            #     os.remove(os.path.join(results, f))

            if f.endswith('.vtk'):
                os.remove(os.path.join(results, f))


    if dd['deletar_npyfiles_flying'] == True:

        flying = 'flying'
        ff = os.listdir(flying)

        for f in ff:
            if f.endswith('.npy'):
                os.remove(os.path.join(flying, f))
