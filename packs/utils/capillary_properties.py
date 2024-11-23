
class capillaryProperties:

    def convert_capillary_pressure_to():
        doc = "The convert_capillary_pressure_to property."
        def fget(self):
            return self._convert_capillary_pressure_to
        def fset(self, value):
            if value not in capillaryPressureBiphasic.names_pressure:
                raise NameError(f'{value} nao esta listado.\nNomes validos: {capillaryPressureBiphasic.names_pressure}')
            self._convert_capillary_pressure_to = value
        def fdel(self):
            # del self._convert_capillary_pressure_to
            pass
        return locals()
    convert_capillary_pressure_to = property(**convert_capillary_pressure_to())
