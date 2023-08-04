from __future__ import print_function, division

import numpy as np
import math
from scipy.stats import truncnorm

import xpsi
from xpsi.global_imports import _G, _csq, _km, _M_s, _2pi
from xpsi.global_imports import gravradius, inv_gravradius

from xpsi.cellmesh.mesh_tools import eval_cedeCentreCoords

from scipy.interpolate import Akima1DInterpolator

a_f = 0.0
b_f = 2.0
a_xi = 0.001
b_xi = math.pi/2.0 - a_xi

class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.

    Currently tailored to the NICER light-curve SWG model specification.

    Source: PSR J0030+0451
    Model variant: ST+PST

    Parameter vector:

    * p[0] = distance (kpc)
    * p[1] = (rotationally deformed) gravitational mass (solar masses)
    * p[2] = coordinate equatorial radius (km)
    * p[3] = inclination of Earth to rotational axis (radians)
    * p[4] = primary centre colatitude (radians)
    * p[5] = primary angular radius (radians)
    * p[6] = primary log10(comoving NSX FIH effective temperature [K])
    * p[7] = secondary centre colatitude (radians)
    * p[8] = secondary angular radius (radians)
    * p[9] = secondary hole colatitude (radians)
    * p[10] = secondary hole angular radius (radians)
    * p[11] = secondary hole azimuth (radians); periodic
    * p[12] = secondary log10(comoving NSX FIH effective temperature [K])
    * p[13] = hydrogen column density (10^20 cm^-2)
    * p[14] = instrument parameter a
    * p[15] = instrument parameter b
    * p[16] = instrument parameter c
    * p[17] = primary cap phase shift (cycles); (alias for initial azimuth, periodic)
    * p[18] = secondary cap phase shift (cycles)

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

        vals = np.linspace(0.0, b_xi, 1000)

        self._interpolator = Akima1DInterpolator(self._vector_super_radius_mass(vals), vals)
        self._interpolator.extrapolate = True

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

        # spots cannot overlap
        theta_p = p[4]
        phi_s = (0.5 + p[18]) * _2pi - p[11]
        phi = p[17] * _2pi - phi_s # include ceding azimuth
        rho_p = p[5]

        theta_s = p[7]
        rho_s = p[8]

        ang_sep = xpsi.Spot._psi(theta_s, phi, theta_p)

        if ang_sep < rho_p + rho_s:
            return -np.inf

        return 0.0

    @staticmethod
    def _I(x):
        return x * np.log(b_xi/a_xi)

    @staticmethod
    def _II(x):
        return 2.0*(x - a_xi) - x*np.log(x/b_xi)

    def _scalar_super_radius_mass(self, x):
        if x >= a_xi:
            mass = self._II(x)
        else:
            mass = self._I(x)

        return mass

    def _vector_super_radius_mass(self, x):
        masses = np.zeros(len(x))

        for i, _ in enumerate(x):
            masses[i] = self._scalar_super_radius_mass(_)

        masses /= (b_f - a_f)
        masses /= (b_xi - a_xi)

        return masses

    @staticmethod
    def _inverse_sample_cede_radius(x, psi):
        if psi < a_xi:
            return a_xi*np.exp(x * np.log(b_xi/a_xi))
        elif psi >= a_xi and x <= 1.0/(1.0 + np.log(b_xi/psi)):
            return x*psi*(1.0 + np.log(b_xi/psi))
        else:
            return psi*np.exp(x*(1.0 + np.log(b_xi/psi)) - 1.0)

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

        # hole radius
        p[10] = float(self._interpolator(hypercube[10]))

        # cede radius
        p[8] = self._inverse_sample_cede_radius(hypercube[8], p[10])

        if p[10] <= p[8]:
            p[7] = hypercube[7] * (p[8] + p[10])
        else:
            p[7] = p[10] - p[8] + 2.0*hypercube[7]*p[8]

        p[7], p[11] = eval_cedeCentreCoords(p[9], p[7], p[11])

        p[11] *= -1.0

        if p[-2] > 0.5:
            p[-2] -= 1.0

        if p[-1] > 0.5:
            p[-1] -= 1.0

        return p

    def inverse_sample_and_transform(self, hypercube):
        """ A transformation for post-processing. """

        p = self.transform(self.inverse_sample(hypercube))

        return p

    @staticmethod
    def transform(p):
        """ A transformation for post-processing. """

        if not isinstance(p, list):
            p = list(p)

        p += [gravradius(p[1]) / p[2]]

        p += [p[8] - p[10]]

        if p[18] > 0.0:
            p += [p[18] - 1.0]
        else:
            p += [p[18]]

        temp = eval_cedeCentreCoords(-1.0*p[9], p[7], -1.0*p[11])

        azi = temp[1]

        if azi < 0.0:
            azi += 2.0*math.pi

        p += [p[10]/p[8] if p[10] <= p[8] else 2.0 - p[8]/p[10]] # f

        p += [p[8] if p[10] <= p[8] else p[10]] # xi

        p += [temp[0]/(p[8] + p[10]) if p[10] <= p[8] else (temp[0] - p[10] + p[8])/(2.0*p[8])] # kappa

        p += [azi/math.pi]

        return p
