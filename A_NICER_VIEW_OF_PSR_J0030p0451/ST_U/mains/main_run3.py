from __future__ import print_function, division

import numpy as np
import math

import xpsi

print('Rank reporting: %d' % xpsi._rank)

from CustomData import CustomData
from CustomInstrument import CustomInstrument
from CustomInterstellar import CustomInterstellar
from CustomPulse import CustomPulse
from CustomSpacetime import CustomSpacetime
from CustomPrior import CustomPrior
from CustomPhotosphere import CustomPhotosphere

data = CustomData.from_SWG('data/NICER_J0030_PaulRay_fixed_evt_25to299__preprocessed.txt', 1936864.0)

NICER = CustomInstrument.from_SWG(num_params=3,
                    bounds=[(0.5,1.5),(0.0,1.0),(0.5,1.5)],
                    ARF = 'model_data/ni_xrcall_onaxis_v1.02_arf.txt',
                    RMF = 'model_data/nicer_upd_d49_matrix.txt',
                    ratio = 'model_data/crab_ratio_SA80_d49.txt',
                    max_input=700,
                    min_input=0,
                    chan_edges = 'model_data/nicer_upd_energy_bounds.txt')

interstellar = CustomInterstellar.from_SWG('model_data/interstellar_phot_frac.txt',
                                           num_params = 1,
                                           bounds = [(0.0, 5.0)])

pulse = CustomPulse(tag = 'all',
                    num_params = 2,
                    bounds = [(-0.25, 0.75), (-0.25, 0.75)],
                    data = data,
                    instrument = NICER,
                    interstellar = interstellar,
                    energies_per_interval = 0.2,
                    workspace_intervals = 1000,
                    epsrel = 1.0e-8,
                    epsilon = 1.0e-3,
                    sigmas = 10.0)

from xpsi.global_imports import _c, _G, _M_s, _dpr, gravradius

bounds = [(0.235, 0.415),
          (1.0, 3.0),
          (3.0 * gravradius(1.0), 16.0),
          (0.001, math.pi/2.0)]

spacetime = CustomSpacetime(num_params = 4, bounds = bounds, S = 1.0/(4.87e-3))

bounds = [(0.001, math.pi - 0.001),
          (0.001, math.pi/2.0 - 0.001),
          (5.1, 6.8),
          (0.001, math.pi - 0.001),
          (0.001, math.pi/2.0 - 0.001),
          (5.1, 6.8)]

spot = xpsi.Spots(num_params=(3,3), bounds=bounds,
                  symmetry=True,
                  cede=False,
                  mirror=False,
                  antipodal=False,
                  translated=False,
                  sq_num_cells=32,
                  num_leaves=64,
                  num_rays = 200)

photosphere = CustomPhotosphere(num_params = 0, bounds = [],
                                tag = 'all', spot = spot, elsewhere = None)

photosphere.spot_atmosphere = 'model_data/nsx_H_v171019.out'

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

likelihood = xpsi.Likelihood(star = star, pulses = pulse, threads=1)

prior = CustomPrior(bounds=likelihood.bounds, spacetime=spacetime)

likelihood.prior = prior

import time

p = [0.338580728497776928E+00,
     0.126592790363188357E+01,
     0.116811279741957499E+02,
     0.103333849972806169E+01,
     0.247170209694014265E+01,
     0.157254430244441135E+00,
     0.607352805901467452E+01,
     0.275378514783974637E+01,
     0.302695201305424211E+00,
     0.606835900698871100E+01,
     0.163598857321082236E+01,
     0.883426307268285194E+00,
     0.329699577808043764E+00,
     0.760519724577322864E+00,
     0.461159408954855521E+00,
     -0.499999E+00]

t = time.time()
ll = likelihood(p) # OptiPlex: ll = -36437.114489180414
print('p: ', ll, time.time() - t)

runtime_params = {'resume': False,
                  'importance_nested_sampling': False,
                  'multimodal': True,
                  'n_clustering_params': None,
                  'outputfiles_basename': './run3_nlive1000_eff0.3_noCONST_MM_tol-1',
                  'n_iter_before_update': 100,
                  'n_live_points': 1000,
                  'sampling_efficiency': 0.3,
                  'const_efficiency_mode': False,
                  'wrapped_params': [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1],
                  'evidence_tolerance': 0.1,
                  'max_iter': -1,
                  'verbose': True}

xpsi.Sample.MultiNest(likelihood, prior, **runtime_params)
