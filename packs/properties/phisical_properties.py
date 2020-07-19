
class PhisicalProperties:

    def __init__(self, input_file=''):
        if input_file == '':
            input_file = 'input_cards/physic/physic_card.yml'

        with open(input_file, 'r') as f:
            data_loaded = yaml.safe_load(f)

        self._gravity = data_loaded['gravity']
        self._gravity_vector = np.array(data_loaded['gravity_vector'], dtype=float)

    @property
    def gravity_vector(self):
        return self._gravity_vector
