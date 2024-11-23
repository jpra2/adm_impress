import numpy as np


def get_req_vec(k00, k11, h0, h1):

    d1 = k00 / k11
    d2 = k11 / k00
    h0_ = h0 * h0
    h1_ = h1 * h1

    temp1 = 0.14 * np.sqrt(np.sqrt(d1) * h0_ + np.sqrt(d2) * h1_)
    temp2 = 0.5 * (np.power(d1, 1 / 4) + np.power(d2, 1 / 4))
    return temp1 / temp2

def get_WI_vec(h, k00, k11, req, well_radius, s=0):
    """

    :param h:
    :param k00: permeability in transversal direction
    :param k11: another permeability in transversal direction
    :param req: equivalent radius
    :param well_radius:
    :return: wi
    """

    wi = 2 * np.pi * ((h * np.sqrt(k00 * k11)) / (np.log(req / np.repeat(well_radius, len(req))) + s))
    return wi
