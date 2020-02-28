import numpy as np


class delta_time:
    def __init__(self, fprop):
        self.P = fprop.P
        self.So = fprop.So
        self.Sg = fprop.Sg
        self.Sw = fprop.Sw
        self.component_mole_numbers = fprop.component_mole_numbers

    def update_deltaT(self, deltaT, fprop):
        """ the limit parameters would be given as data entry """
        deltaPlim = 1
        deltaSlim = 0.1
        deltaNlim = 0.01
        deltaVlim = 0.01

        ''' Still confused of how I will calculate this maximum deltas if they
        depend of the new time step, who depends of deltaT...'''

        deltaPmax = max(np.abs(fprop.P - self.P) / fprop.P)
        deltaSomax = np.nanmax(np.abs(fprop.So - self.So))
        deltaSgmax = np.nanmax(np.abs(fprop.Sg - self.Sg))
        deltaSwmax = np.nanmax(np.abs(fprop.Sw - self.Sw))
        deltaSmax = max(deltaSomax, deltaSgmax, deltaSwmax)

        old_settings = np.seterr(all = 'ignore', divide = 'ignore')
        deltaNmax = max(np.nanmax(np.abs(fprop.component_mole_numbers - self.component_mole_numbers)
                        / fprop.component_mole_numbers, axis =1))

        np.seterr(**old_settings)
        deltaVmax = max(np.abs(fprop.Vt - fprop.Vp) / fprop.Vp)

        deltaTp = deltaT * deltaPlim / deltaPmax
        deltaTs = deltaT * deltaSlim / deltaSmax
        deltaTn = deltaT * deltaNlim / deltaNmax
        deltaTv = deltaT * deltaVlim / deltaVmax
        deltaT = min(deltaTp, deltaTs, deltaTn, deltaTv)
        return deltaT
