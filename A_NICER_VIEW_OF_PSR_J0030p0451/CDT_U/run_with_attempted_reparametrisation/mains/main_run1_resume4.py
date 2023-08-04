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
                    energies_per_interval = 0.25,
                    fast_rel_energies_per_interval = 0.5,
                    default_energy_spacing = 'logspace',
                    adaptive_energies = False,
                    adapt_exponent = None,
                    store = False,
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
          (0.001, math.pi/2.0 - 0.001),
          (5.1, 6.8),
          (5.1, 6.8),
          (0.001, math.pi - 0.001),
          (0.001, math.pi/2.0 - 0.001),
          (0.001, math.pi/2.0 - 0.001),
          (5.1, 6.8),
          (5.1, 6.8)]

spot = xpsi.Spots(num_params=(5,5), bounds=bounds,
                  symmetry=False,
                  cede=True,
                  concentric=True,
                  antipodal_symmetry=False,
                  sqrt_num_cells=24,
                  min_sqrt_num_cells=10,
                  num_leaves=80,
                  num_rays=200,
                  fast_sqrt_num_cells=8,
                  fast_min_sqrt_num_cells=8,
                  fast_num_leaves=32,
                  fast_num_rays=100)

photosphere = CustomPhotosphere(num_params = 0, bounds = [],
                                tag = 'all', spot = spot, elsewhere = None)

photosphere.spot_atmosphere = 'model_data/nsx_H_v171019.out'

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

likelihood = xpsi.Likelihood(star = star, pulses = pulse,
                             fast_precomputation = True, threads=1)

prior = CustomPrior(bounds=likelihood.bounds, spacetime=spacetime)

likelihood.prior = prior

import time
p = [0.326524654427988337E+00,
        0.148829365719602813E+01,
        0.137747061119835159E+02,
        0.933292525886352808E+00,
        0.217467609327838263E+01,
        0.267191210890304620E+00 * 0.261151597963538817E+00,
        0.267191210890304620E+00,
        0.610946708456508780E+01,
        0.553945149821494187E+01,
        0.251371132748997850E+01,
        0.310173011684981348E+00 * 0.927951102138210637E+00,
        0.310173011684981348E+00,
        0.570882469885388577E+01,
        0.612226554815461466E+01,
        0.106317116410402046E+01,
        0.110313346338052098E+01,
        0.479840116654693463E-01,
        0.115782124307722700E+01,
        0.454366838390029781E+00,
        0.496943822277659319E+00]

t = time.time()
ll = likelihood(p) # OptiPlex: ll = -36317.23819085545
print('p: ', ll, time.time() - t)

runtime_params = {'resume': True,
                  'importance_nested_sampling': False,
                  'multimodal': False,
                  'n_clustering_params': None,
                  'outputfiles_basename': './run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1',
                  'n_iter_before_update': 50,
                  'n_live_points': 1000,
                  'sampling_efficiency': 0.8,
                  'const_efficiency_mode': False,
                  'wrapped_params': [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1],
                  'evidence_tolerance': 0.1,
                  'max_iter': -1,
                  'verbose': True}

xpsi.Sample.MultiNest(likelihood, prior, **runtime_params)
