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
                    num_params = 1,
                    bounds = [(-0.25, 0.75)],
                    data = data,
                    instrument = NICER,
                    interstellar = interstellar,
                    energies_per_interval = 0.25,
                    default_energy_spacing = 'logspace',
                    fast_rel_energies_per_interval = 0.5,
                    workspace_intervals = 1000,
                    adaptive_energies=False,
                    adapt_exponent=0.5,
                    store=False,
                    epsrel = 1.0e-8,
                    epsilon = 1.0e-3,
                    sigmas = 10.0)

from xpsi.global_imports import _c, _G, _M_s, _dpr, gravradius

bounds = [(0.235, 0.415),
          (1.0, 3.0),
          (3.0 * gravradius(1.0), 16.0),
          (0.001, math.pi/2.0)]

spacetime = CustomSpacetime(num_params = 4, bounds = bounds, S = 1.0/(4.87e-3))

bounds = [(0.001, math.pi/2.0),
          (0.001, math.pi/2.0 - 0.001),
          (5.1, 6.8)]

spot = xpsi.Spots(num_params=(3,0), bounds=bounds,
                    symmetry=True,
                    cede=False,
                    hole=False,
                    concentric=False,
                    antipodal_symmetry=True,
                    sqrt_num_cells=32,
                    min_sqrt_num_cells=10,
                    max_sqrt_num_cells=64,
                    do_fast=False,
                    fast_sqrt_num_cells=8,
                    fast_min_sqrt_num_cells=8,
                    fast_max_sqrt_num_cells=16,
                    fast_num_leaves=32,
                    fast_num_rays=100,
                    num_leaves=100,
                    num_rays=200,
                    is_secondary=False)

photosphere = CustomPhotosphere(num_params = 0, bounds = [],
                                tag = 'all', spot = spot, elsewhere = None)

photosphere.spot_atmosphere = 'model_data/nsx_H_v171019.out'

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

likelihood = xpsi.Likelihood(star = star, pulses = pulse, threads=1)

prior = CustomPrior(bounds=likelihood.bounds, spacetime=spacetime)

likelihood.prior = prior

import time

p = [0.319405461347693653E+00,
        0.293767931853114073E+01,
        0.159997443154808501E+02,
        0.122161682172283914E+01,
        0.126378196054516656E+01,
        0.101144501416911836E+00,
        0.604435138549768958E+01,
        0.595127021740757527E-02,
        0.104664618877441340E+01,
        0.161578048557822540E+00,
        0.997221406820496759E+00,
        -0.882459875267529448E-01]

t = time.time()
ll = likelihood(p) # OptiPlex: ll = -37261.546911809055
print('p: ', ll, time.time() - t)

runtime_params = {'resume': False,
                  'importance_nested_sampling': False,
                  'multimodal': False,
                  'n_clustering_params': None,
                  'outputfiles_basename': './run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1',
                  'n_iter_before_update': 100,
                  'n_live_points': 1000,
                  'sampling_efficiency': 0.3,
                  'const_efficiency_mode': False,
                  'wrapped_params': [0,0,0,0,0,0,0,0,0,0,0,1],
                  'evidence_tolerance': 0.1,
                  'max_iter': -1,
                  'verbose': True}

xpsi.Sample.MultiNest(likelihood, prior, **runtime_params)
