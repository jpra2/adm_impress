import numpy as np


class delta_time:
    def __init__(self, fprop):
        self.P = fprop.P
        self.So = fprop.So
        self.Sg = fprop.Sg
        self.Sw = fprop.Sw
        self.component_mole_numbers = fprop.component_mole_numbers
        #the initialization of this class is made in a different time step evaluation

    def update_CFL(deltaT, fprop):
        old_settings = np.seterr(all = 'ignore', divide = 'ignore')
        CFL = np.nanmax(deltaT * fprop.component_flux_vols_total[fprop.component_mole_numbers!=0] /
        fprop.component_mole_numbers[fprop.component_mole_numbers!=0])
        if (CFL > 1): deltaT = deltaT / 2

        #deltaTcfl = np.nanmin(CFL * (fprop.component_mole_numbers / fprop.Vbulk) / fprop.component_flux_vols_total, axis = 1) #make nan
        np.seterr(**old_settings)
        return deltaT

    def update_deltaTp(self, deltaT, fprop, deltaPlim):
        deltaPmax = max(np.abs(fprop.P - self.P) / fprop.P)
        deltaTp = deltaT * deltaPlim / deltaPmax
        return deltaTp

    def update_deltaTs(self, deltaT, fprop, deltaSlim):
        old_settings = np.seterr(all = 'ignore', divide = 'ignore')
        deltaSo = np.abs(fprop.So - self.So) / fprop.So
        deltaSg = np.abs(fprop.Sg - self.Sg) / fprop.Sg
        deltaSw = np.abs(fprop.Sw - self.Sw) / fprop.Sw

        deltasS = np.array([deltaSo,deltaSg,deltaSw])
        deltaSmax = np.nanmax(deltasS)
        np.seterr(**old_settings)

        deltaTs = deltaT * deltaSlim / deltaSmax
        return deltaTs

    def update_deltaTn(self, deltaT, fprop, deltaNlim):
        old_settings = np.seterr(all = 'ignore', divide = 'ignore')
        deltaNmax = max(np.nanmax(np.abs(fprop.component_mole_numbers - self.component_mole_numbers)
                        / fprop.component_mole_numbers, axis =1))
        np.seterr(**old_settings)
        deltaTn = deltaT * deltaNlim / deltaNmax
        return deltaTn

    def update_deltaTv(self, deltaT, fprop, deltaVlim):
        deltaVmax = max(np.abs(fprop.Vt - fprop.Vp) / fprop.Vp)
        deltaTv = deltaT * deltaVlim / deltaVmax
        return deltaTv

    def update_deltaT(self, deltaT, fprop):
        """ the limit parameters would be given as data entry """
        deltaPlim = 1
        deltaSlim = 0.1
        deltaNlim = 0.01
        deltaVlim = 0.01

        ''' Still confused of how I will calculate this maximum deltas if they
        depend of the new time step, who depends of deltaT...'''

        deltaTp = self.update_deltaTp(deltaT, fprop, deltaPlim)
        deltaTs = self.update_deltaTs(deltaT, fprop, deltaSlim)
        deltaTn = self.update_deltaTs(deltaT, fprop, deltaNlim)
        deltaTv = self.update_deltaTs(deltaT, fprop, deltaVlim)

        deltaT = min(deltaTp, deltaTs, deltaTn, deltaTv)
        return deltaT
