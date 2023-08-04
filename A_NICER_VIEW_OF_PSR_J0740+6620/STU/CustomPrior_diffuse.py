from __future__ import print_function, division

import numpy as np
import math
from scipy.stats import truncnorm, skewnorm
from scipy.integrate import quad
from scipy.interpolate import Akima1DInterpolator

import xpsi
from xpsi.global_imports import _2pi, gravradius, _dpr
from xpsi import Parameter

_cond_std_alpha  = np.sqrt( 0.75 * 0.15 * 0.15 )
l = 149.7 * math.pi/180.0
b = 29.60 * math.pi/180.0
R_o = 8.5
H = 0.2 * R_o
h = 0.5

def distance_PDF(D):
    global l
    global b

    R = R_o*R_o
    R += np.sqrt(R_o*R_o + (D*math.cos(b))**2.0 - 2.0*D*math.cos(b)*R_o*math.cos(l))
    z = D * np.sin(b)

    F_D = R**1.9 * np.exp(-1.0 * np.fabs(z)/h - R/H)

    #F_D *= np.exp(-0.5 * ((1.0/D - 2.3)/0.6)**2.0)

    return D * D * F_D / 6.811441705752484e-18

distances = np.linspace(0.1, 10.0, 1000)
CDF = np.zeros(1000)
for i in range(len(distances)):
    CDF[i] = quad(distance_PDF, 0.1, distances[i])[0]

interpolator = Akima1DInterpolator(CDF, distances)
interpolator.extrapolate = True

class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.

    Source: PSR J0740
    Model variant: ST-U
        Two single-temperature, simply-connected circular hot regions with
        unshared parameters.

    Parameter vector: (print the likelihood object)

    * p[0] = (rotationally deformed) gravitational mass (solar masses)
    * p[1] = coordinate equatorial radius (km)
    * p[2] = distance (kpc)
    * p[3] = cos(inclination of Earth to rotational axis)
    * p[4] = primary cap phase shift (cycles); (alias for initial azimuth, periodic)
    * p[5] = primary centre colatitude (radians)
    * p[6] = primary angular radius (radians)
    * p[7] = primary log10(comoving NSX FIH effective temperature [K])
    * p[8] = secondary cap phase shift (cycles)
    * p[9] = secondary centre colatitude (radians)
    * p[10] = secondary angular radius (radians)
    * p[11] = secondary log10(comoving NSX FIH effective temperature [K])
    * p[12] = energy-independent NICER XTI effective-area scaling
    * p[13] = hydrogen column density (10^20 cm^-2)
    * p[14] = energy-independent XMM effective-area scaling (relative to XTI)

    """

    __derived_names__ = ['compactness',
                         'p__phase_shift_shifted',
                         's__phase_shift_shifted']

    __draws_from_support__ = 4

    def __call__(self, p = None):
        """ Evaluate distribution at ``p``.

        :param list p: Model parameter values.

        :returns: Logarithm of the distribution evaluated at ``p``.

        """
        temp = super(CustomPrior, self).__call__(p)
        if not np.isfinite(temp):
            return temp

        #if not 0.0 < self.parameters['distance'] <= 1.7: # Shklovskii trunc.
        #    return -np.inf

        # based on contemporary EOS theory
        if not self.parameters['radius'] <= 16.0:
            return -np.inf

        ref = self.parameters.star.spacetime # shortcut

        # limit polar radius to be outside the Schwarzschild photon sphere
        R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
        if R_p < 1.505 / ref.R_r_s:
            return -np.inf

        mu = math.sqrt(-1.0 / (3.0 * ref.epsilon * (-0.788 + 1.030 * ref.zeta)))

        # 2-surface cross-section have a single maximum in |z|
        # i.e., an elliptical surface; minor effect on support, if any,
        # only for high spin frequencies
        if mu < 1.0:
            return -np.inf

        # check effective gravity at pole (where it is maximum) and
        # at equator (where it is minimum) are in NSX limits
        grav = xpsi.surface_radiation_field.effective_gravity(np.array([1.0, 0.0]),
                                                              np.array([ref.R] * 2 ),
                                                              np.array([ref.zeta] * 2),
                                                              np.array([ref.epsilon] * 2))
        for g in grav:
            if not 13.7 <= g <= 15.0:
                return -np.inf

        ref = self.parameters # redefine shortcut

        # enforce order in hot region colatitude
        if ref['p__super_colatitude'] > ref['s__super_colatitude']:
            return -np.inf

        phi = (ref['p__phase_shift'] - 0.5 - ref['s__phase_shift']) * _2pi

        ang_sep = xpsi.HotRegion.psi(ref['s__super_colatitude'],
                                     phi,
                                     ref['p__super_colatitude'])

        # hot regions cannot overlap
        if ang_sep < ref['p__super_radius'] + ref['s__super_radius']:
            return -np.inf

        return 0.0

    def inverse_sample(self, hypercube=None):
        """ Draw sample uniformly from the distribution via inverse sampling. """
        global interpolator

        to_cache = self.parameters.vector

        if hypercube is None:
            hypercube = np.random.rand(len(self))

        # the base method is useful, so to avoid writing that code again:
        _ = super(CustomPrior, self).inverse_sample(hypercube)

        ref = self.parameters # redefine shortcut

        idx = ref.index('distance')
        ref['distance'] = interpolator(hypercube[idx])

        #idx = ref.index('mass')
        #ref['mass'] = truncnorm.ppf(hypercube[idx], -5.0, 5.0,
        #                            loc= , scale= )

        idx = ref.index('cos_inclination')
        ref['cos_inclination'] = hypercube[idx]

        # flat priors in cosine of hot region centre colatitudes (isotropy)
        # support modified by no-overlap rejection condition
        idx = ref.index('p__super_colatitude')
        a, b = ref.get_param('p__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['p__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])

        idx = ref.index('s__super_colatitude')
        a, b = ref.get_param('s__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['s__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])

        try:
            ref['XTI__alpha']
        except KeyError:
            pass
        else:
            idx = ref.index('XTI__alpha')
            ref['XTI__alpha'] = truncnorm.ppf(hypercube[idx], -5.0, 5.0,
                                          loc=1.0, scale=0.15)

        try:
            ref['PN__alpha']
        except KeyError:
            pass
        else:
            idx = ref.index('PN__alpha')
            try:
                ref['XTI__alpha']
            except KeyError:
                _loc = 1.0
                _scale = 0.15
            else:
                _loc = 1.0 + 0.5 * (ref['XTI__alpha'] - 1.0)
                _scale = _cond_std_alpha
            ref['PN__alpha'] = truncnorm.ppf(hypercube[idx], -5.0, 5.0,
                                         loc=_loc,
                                         scale=_scale)

        #ref['absolute_alpha'] = truncnorm.ppf(hypercube[ ref.index('absolute_alpha') ],
        #                             -5.0, 5.0, loc=1.0, scale=0.1)

        #for instrument in ['XTI', 'PN']:
        #    _name = instrument + '__alpha'
        #    try:
        #        ref[_name] = truncnorm.ppf(hypercube[ ref.index(_name) ],
        #                                      -5.0, 5.0, loc=1.0, scale=0.1)
        #    except KeyError:
        #        pass

        # restore proper cache
        for parameter, cache in zip(self.parameters, to_cache):
            parameter.cached = cache

        # it is important that we return the desired vector because it is
        # automatically written to disk by MultiNest and only by MultiNest
        return self.parameters.vector

    def transform(self, p, **kwargs):
        """ A transformation for post-processing. """

        p = list(p) # copy

        # used ordered names and values
        ref = dict(zip(self.parameters.names, p))

        # compactness ratio M/R_eq
        p += [gravradius(ref['mass']) / ref['radius']]

        for phase_shift in ['p__phase_shift', 's__phase_shift']:
            if ref[phase_shift] > 0.5:
                p += [ref[phase_shift] - 1.0]
            else:
                p += [ref[phase_shift]]

        return p
