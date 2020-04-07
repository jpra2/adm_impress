
class biphasicProperties:

    def gama_volumes_average():
        doc = "The gama_volumes_average property."
        def fget(self):
            gama_w = self.biphasic_data['gama_w']
            gama_o = self.biphasic_data['gama_o']
            lambda_w = self.data_impress['lambda_w']
            lambda_o = self.data_impress['lambda_o']
            return np.repeat(gama_w, len(lambda_w))*lambda_w + np.repeat(gama_o, len(lambda_o))*lambda_o
        return locals()
    gama_volumes_average = property(**gama_volumes_average())

    def lambda_w_internal_faces():
        doc = "The lambda_w_internal_faces property."
        def fget(self):
            # return self.data_impress['lambda_w'][self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            return self.lambda_w_volumes[self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
        return locals()
    lambda_w_internal_faces = property(**lambda_w_internal_faces())

    def lambda_o_internal_faces():
        doc = "The lambda_o_internal_faces property."
        def fget(self):
            # return self.data_impress['lambda_o'][self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            return self.lambda_o_volumes[self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
        return locals()
    lambda_o_internal_faces = property(**lambda_o_internal_faces())

    def lambda_t_internal_faces():
        doc = "The lambda_t_internal_faces property."
        def fget(self):
            # return self.data_impress['lambda_t'][self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            return self.lambda_t_volumes[self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
        return locals()
    lambda_t_internal_faces = property(**lambda_t_internal_faces())

    def fw_internal_faces():
        doc = "The fw_internal_faces property."
        def fget(self):
            # return self.data_impress['fw_vol'][self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            return self.fw_volumes[self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
        return locals()
    fw_internal_faces = property(**fw_internal_faces())

    def lambda_w_volumes():
        doc = "The lambda_w_volumes property."
        def fget(self):
            return self.data_impress['krw']/self.biphasic_data['mi_w']
        return locals()
    lambda_w_volumes = property(**lambda_w_volumes())

    def lambda_o_volumes():
        doc = "The lambda_o_volumes property."
        def fget(self):
            return self.data_impress['kro']/self.biphasic_data['mi_o']
        return locals()
    lambda_o_volumes = property(**lambda_o_volumes())

    def lambda_t_volumes():
        doc = "The lambda_t_volumes property."
        def fget(self):
            return self.lambda_w_volumes + self.lambda_o_volumes
        return locals()
    lambda_t_volumes = property(**lambda_t_volumes())

    def fw_volumes():
        doc = "The fw_volumes property."
        def fget(self):
            return self.lambda_w_volumes/self.lambda_t_volumes
        return locals()
    fw_volumes = property(**fw_volumes())
