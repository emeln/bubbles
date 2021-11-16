from .sigma_interpolation import SigmaInterpolator
from .cosmology import Cosmology
from scipy.special import erfinv
import numpy as np


class BubbleDist(object):

    def __init__(self, cosmo, sigmaInterpolator):
        """

                :type sigmaInterpolator: SigmaInterpolator
                :type cosmo: Cosmology
                """
        self.cosmo = cosmo
        self.sigmaInt = sigmaInterpolator
        self.__zeta = 40


    def __call__(self, *args, **kwargs):
        raise NotImplementedError("This is not implemented here! Use sublcasses!")

    @property
    def zeta(self):
        return self.__zeta

    @zeta.Setter
    def zeta(self, new):
        self.__zeta = new


class BubbleDistributionFunction(BubbleDist):
    def __call__(self, m, z):
        """

        output is dN/d log m (= m dN/dm ??)

        In Furlanetto, delta_crit is a function of z, but Garrit (& Bartelmann?) have it as a
        constant delta_crit ~ 1.69
        """
        mMin = 1.38e6 * (kb / (G * mp))**(3/2) * (overdensity(z)*)

        vals = np.sqrt(2 / np.pi) * (self.cosmo.rho_mean / m ) * np.fabs(
                self.sigmaInt.dlogSigma_dlogm(m, z)) * (self.__B0(z) / self.sigmaInt(m, z)) * np.exp(-self.__B(m, z) ** 2 / (2 * self.sigmaInt(m, z) ** 2))

        # cosmo.rho_mean returns mean mass density today -- should that have z dependence in this equation?

        return vals


    def __B0(self, z):
        return self.cosmo.delta_crit - np.sqrt(2) * erfinv(1 - 1/self.__zeta) * self.sigmaInt(self.cosmo.mVir(z), z)

    def __B(self, m, z):
        return self.__B0(z) + erfinv(1 - 1/self.__zeta) / np.sqrt(2 * (self.sigmaInt(self.cosmo.mVir(z),z)**2 - self.sigmaInt(m,z)**2)) * self.sigmaInt(m,z)**2



