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
                    bounds = [(0.35, 0.55), (-0.25,0.75)],
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
          (5.1, 6.8)]

primary = xpsi.Spot(num_params=3, bounds=bounds,
                    symmetry=True,
                    hole=False,
                    cede=False,
                    concentric=False,
                    sqrt_num_cells=24,
                    min_sqrt_num_cells=10,
                    max_sqrt_num_cells=64,
                    do_fast=False,
                    fast_sqrt_num_cells=8,
                    fast_min_sqrt_num_cells=8,
                    fast_max_sqrt_num_cells=16,
                    fast_num_leaves=32,
                    fast_num_rays=100,
                    num_leaves=80,
                    num_rays=200)

bounds = [(0.001, math.pi - 0.001),
          (0.001, math.pi/2.0 - 0.001),
          (0.001, math.pi - 0.001),
          (0.0, 1.0),
          (0.0, 2.0*math.pi),
          (5.1, 6.8)]

secondary = xpsi.Spot(num_params=6, bounds=bounds,
                      symmetry=True,
                      hole=True,
                      cede=False,
                      concentric=False,
                      sqrt_num_cells=24,
                      min_sqrt_num_cells=10,
                      max_sqrt_num_cells=64,
                      do_fast=False,
                      fast_sqrt_num_cells=8,
                      fast_min_sqrt_num_cells=8,
                      fast_max_sqrt_num_cells=16,
                      fast_num_leaves=32,
                      fast_num_rays=100,
                      num_leaves=80,
                      num_rays=200,
                      is_secondary=True)

from xpsi import TwoSpots

spot = TwoSpots((primary, secondary))

photosphere = CustomPhotosphere(num_params = 0, bounds = [],
                                tag = 'all', spot = spot, elsewhere = None)

photosphere.spot_atmosphere = 'model_data/nsx_H_v171019.out'

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

likelihood = xpsi.Likelihood(star = star, pulses = pulse, threads=1)

prior = CustomPrior(bounds=likelihood.bounds, spacetime=spacetime)

likelihood.prior = prior

import time

p = [0.336793253850493635E+00,
        0.159242252411888296E+01,
        0.149436612926094980E+02,
        0.984343441146659059E+00,
        0.223158695736117929E+01,
        0.718883403004594718E-01,
        0.610860748136519405E+01,
        0.258818649563776981E+01,
        0.272499290209381573E+00,
        0.258818649563776981E+01,
        0.237096325731429775E+00,
        0.0,
        0.611204195144800799E+01,
        0.613859545134022522E+00,
        0.971055676486793806E+00,
        0.213885526029857452E-03,
        0.994127256659159242E+00,
        0.454115827418820062E+00,
        0.498373231219062185E+00]

t = time.time()
ll = likelihood(p) # OptiPlex: ll = -36318.56177006986
print('p: ', ll, time.time() - t)

runtime_params = {'resume': True,
                  'importance_nested_sampling': False,
                  'multimodal': False,
                  'n_clustering_params': None,
                  'outputfiles_basename': './run3_nlive1000_eff0.3_noCONST_noMM_IS_tol-1',
                  'n_iter_before_update': 100,
                  'n_live_points': 1000,
                  'sampling_efficiency': 0.3,
                  'const_efficiency_mode': False,
                  'wrapped_params': [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1],
                  'evidence_tolerance': 0.1,
                  'max_iter': -1,
                  'verbose': True}

xpsi.Sample.MultiNest(likelihood, prior, **runtime_params)
