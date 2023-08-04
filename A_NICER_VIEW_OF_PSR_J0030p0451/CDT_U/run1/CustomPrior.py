from __future__ import print_function, division

import numpy as np
import math
from scipy.stats import truncnorm

import xpsi
from xpsi.global_imports import _G, _csq, _km, _M_s, _2pi

class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.

    Currently tailored to the NICER light-curve SWG model specification.

    Source: PSR J0030+0451
    Model variant: CDT-U

    Parameter vector:

    * p[0] = distance (kpc)
    * p[1] = (rotationally deformed) gravitational mass (solar masses)
    * p[2] = coordinate equatorial radius (km)
    * p[3] = inclination of Earth to rotational axis (radians)
    * p[4] = primary-cede centre colatitude (radians)
    * p[5] = primary-cede angular radius (radians)
    * p[6] = primary-super fractional angular radius
    * p[7] = primary-cede log10(comoving NSX FIH effective temperature [K])
    * p[8] = primary-super log10(comoving NSX FIH effective temperature [K])
    * p[9] = secondary-cede centre colatitude (radians)
    * p[10] = secondary-cede angular radius (radians)
    * p[11] = secondary-super fractional angular radius
    * p[12] = secondary-cede log10(comoving NSX FIH effective temperature [K])
    * p[13] = secondary-super log10(comoving NSX FIH effective temperature [K])
    * p[14] = hydrogen column density (10^20 cm^-2)
    * p[15] = instrument parameter a
    * p[16] = instrument parameter b
    * p[17] = instrument parameter c
    * p[18] = primary cap phase shift (cycles); (alias for initial azimuth, periodic)
    * p[19] = secondary cap phase shift (cycles)

    Note that the unit hypercube to physical transformation is constructed
    for the phases by inverse sampling a flat prior on [-0.25,0.75].
    There is then no need for a periodic boundary and we need to worry about
    accuracy at the boundary.

    """
    def __init__(self, bounds, spacetime):
        # Execute abstract parent initialiser
        super(CustomPrior, self).__init__(bounds)

        assert isinstance(spacetime, xpsi.Spacetime),\
                'Invalid type for ambient spacetime object.'

        self._spacetime = spacetime

    def __call__(self, p):
        """ Evaluate distribution at :obj:`p`.

        :param list p: Model parameters values.

        :return: Logarithm of the distribution evaluated at :obj:`p`.

        """
        i = self._spacetime.num_params
        self._spacetime.update(*p[:i])

        if not self._spacetime.R <= 16.0*_km:
            return -np.inf

        if not 1.5 < self._spacetime.R_r_s:
            return -np.inf

        epsilon = self._spacetime.epsilon
        zeta = self._spacetime.zeta
        mu = math.sqrt(-1.0 / (3.0 * epsilon * (-0.788 + 1.030 * zeta)))

        # 2-surface cross-section have a single maximum in |z|
        # i.e., an elliptical surface
        if mu < 1.0:
            return -np.inf

        # polar radius causality for ~static star (static ambient spacetime)
        R_p = 1.0 + epsilon * (-0.788 + 1.030 * zeta)

        if R_p < 1.5 / self._spacetime.R_r_s:
            return -np.inf

        if p[4] > p[9]:
            return -np.inf

        # spots cannot overlap
        theta_p = p[4]
        phi = (p[-2] - 0.5 - p[-1]) * _2pi
        rho_p = p[5]

        theta_s = p[9]
        rho_s = p[10]

        ang_sep = xpsi.Spot._psi(theta_s, phi, theta_p)

        if ang_sep < rho_p + rho_s:
            return -np.inf

        return 0.0

    def inverse_sample(self, hypercube):
        """ Draw sample uniformly from the distribution via inverse sampling.

        :param hypercube: A pseudorandom point in an n-dimensional hypercube.

        :return: A parameter ``list``.

        """
        p = super(CustomPrior, self).inverse_sample(hypercube)

        # distance
        p[0] = truncnorm.ppf(hypercube[0], -10.0, 10.0, loc=0.325, scale=0.009)

        # instrument parameter a
        p[-5] = truncnorm.ppf(hypercube[-5], -5.0, 5.0, loc=1.0, scale=0.1)

        # instrument parameter c
        p[-3] = truncnorm.ppf(hypercube[-3], -5.0, 5.0, loc=1.0, scale=0.1)

        if p[-2] > 0.5:
            p[-2] -= 1.0

        if p[-1] > 0.5:
            p[-1] -= 1.0

        return p
