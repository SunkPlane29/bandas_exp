""" Main module for NICER J0740 <- X-PSI v0.7 ST-U. """
from __future__ import print_function, division

import argparse

parser = argparse.ArgumentParser(
    description='''
    Main module for X-PSI ST-U modelling of NICER J0740+6620 event data.

    You can run this module as a script and launch a sampler, optionally
    with a world of MPI processes.

    Alternate usage: mpiexec -n 4 python -m mpi4py %(prog)s [-h] @<config.ini> [--multinest] [--emcee]

    ''',
    fromfile_prefix_chars='@')

_prfx = 'Absolute or relative path to '
_instr = 'NICER-XTI '
parser.add_argument('--NICER-matrix-path', type=str, help=_prfx + _instr + 'channel-phase count matrix. This path is written to if the file does not exist.')
parser.add_argument('--NICER-event-path', type=str, help=_prfx + _instr + 'event list file.')
parser.add_argument('--NICER-arf-path', type=str, help=_prfx + _instr + 'ARF file.')
parser.add_argument('--NICER-rmf-path', type=str, help=_prfx + _instr + 'RMF file.')
parser.add_argument('--NICER-channels-path', type=str, help=_prfx + _instr + 'channel bounds file.')

_instr = 'XMM-PN '
parser.add_argument('--PN-spectrum-path', type=str, help=_prfx + _instr + 'spectrum file.')
parser.add_argument('--PN-arf-path', type=str, help=_prfx + _instr + 'ARF file.')
parser.add_argument('--PN-rmf-path', type=str, help=_prfx + _instr + 'RMF file.')
parser.add_argument('--PN-channels-path', type=str, help=_prfx + _instr + 'channel bounds file.')
parser.add_argument('--PN-background-path', type=str, help=_prfx + _instr + 'background spectrum file.')

_instr = 'XMM-MOS1 '
parser.add_argument('--MOS1-spectrum-path', type=str, help=_prfx + _instr + 'spectrum file.')
parser.add_argument('--MOS1-arf-path', type=str, help=_prfx + _instr + 'ARF file.')
parser.add_argument('--MOS1-rmf-path', type=str, help=_prfx + _instr + 'RMF file.')
parser.add_argument('--MOS1-channels-path', type=str, help=_prfx + _instr + 'channel bounds file.')
parser.add_argument('--MOS1-background-path', type=str, help=_prfx + _instr + 'background spectrum file.')

_instr = 'XMM-MOS2 '
parser.add_argument('--MOS2-spectrum-path', type=str, help=_prfx + _instr + 'spectrum file.')
parser.add_argument('--MOS2-arf-path', type=str, help=_prfx + _instr + 'ARF file.')
parser.add_argument('--MOS2-rmf-path', type=str, help=_prfx + _instr + 'RMF file.')
parser.add_argument('--MOS2-channels-path', type=str, help=_prfx + _instr + 'channel bounds file.')
parser.add_argument('--MOS2-background-path', type=str, help=_prfx + _instr + 'background spectrum file.')

parser.add_argument('--attenuation-path', type=str, help=_prfx + 'attenuation file.')
parser.add_argument('--atmosphere-path', type=str, help=_prfx + 'atmosphere file.')

parser.add_argument('--multinest', action='store_true',
                    help='Launch MultiNest sampler. Takes precedence.')
parser.add_argument('--emcee', action='store_true',
                    help='Launch emcee sampler.')
parser.add_argument('--resume', action='store_true',
                    help='Resume sampling?')

parser.add_argument('--NICER', action='store_true',
                    help='Model NICER event data.')
parser.add_argument('--XMM', action='store_true',
                    help='Model XMM event data.')

if __name__ == '__main__':
    args = parser.parse_args()
else:
    args = parser.parse_args(['@STU/config_He.ini'])

# check if interactive input needed
# NICER
if not args.NICER_matrix_path:
    args.NICER_matrix_path = raw_input('Specify the NICER response matrix path: ')

if not args.NICER_event_path:
    args.NICER_event_path = raw_input('Specify the NICER event file path: ')

if not args.NICER_arf_path:
    args.NICER_arf_path = raw_input('Specify the NICER ARF file path: ')

if not args.NICER_rmf_path:
    args.NICER_rmf_path = raw_input('Specify the NICER RMF file path: ')

if not args.NICER_channels_path:
    args.NICER_channels_path = raw_input('Specify the NICER channel energy file path: ')

# PN
if not args.PN_spectrum_path:
    args.PN_spectrum_path = raw_input('Specify the XMM-PN spectrum file path: ')

if not args.PN_arf_path:
    args.PN_arf_path = raw_input('Specify the XMM-PN ARF file path: ')

if not args.PN_rmf_path:
    args.PN_rmf_path = raw_input('Specify the XMM-PN RMF file path: ')

if not args.PN_channels_path:
    args.PN_channels_path = raw_input('Specify the XMM-PN channel energy file path: ')

if not args.PN_background_path:
    args.PN_background_path = raw_input('Specify the XMM-PN background file path: ')

# MOS1
if not args.MOS1_spectrum_path:
    args.MOS1_spectrum_path = raw_input('Specify the XMM-MOS1 spectrum file path: ')

if not args.MOS1_arf_path:
    args.MOS1_arf_path = raw_input('Specify the XMM-MOS1 ARF file path: ')

if not args.MOS1_rmf_path:
    args.MOS1_rmf_path = raw_input('Specify the XMM-MOS1 RMF file path: ')

if not args.MOS1_channels_path:
    args.MOS1_channels_path = raw_input('Specify the XMM-MOS1 channel energy file path: ')

if not args.MOS1_background_path:
    args.MOS1_background_path = raw_input('Specify the XMM-MOS1 background file path: ')

# MOS2
if not args.MOS2_spectrum_path:
    args.MOS2_spectrum_path = raw_input('Specify the XMM-MOS2 spectrum file path: ')

if not args.MOS2_arf_path:
    args.MOS2_arf_path = raw_input('Specify the XMM-MOS2 ARF file path: ')

if not args.MOS2_rmf_path:
    args.MOS2_rmf_path = raw_input('Specify the XMM-MOS2 RMF file path: ')

if not args.MOS2_channels_path:
    args.MOS2_channels_path = raw_input('Specify the XMM-MOS2 channel energy file path: ')

if not args.MOS2_background_path:
    args.MOS2_background_path = raw_input('Specify the XMM-MOS2 background file path: ')

# shared
if not args.attenuation_path:
    args.attenuation_path = raw_input('Specify the attenuation file path: ')

if not args.atmosphere_path:
    args.atmosphere_path = raw_input('Specify the atmosphere file path: ')

import numpy as np
import math

import xpsi
from xpsi.Parameter import Derive

#print('Rank reporting: %d' % xpsi._rank)

from xpsi.global_imports import gravradius

from STU.CustomInstrument import CustomInstrument
from STU.CustomSignal import CustomSignal
from STU.CustomInterstellar import CustomInterstellar
from STU.CustomPhotosphere_He import CustomPhotosphere
from STU.CustomPrior import CustomPrior

class namespace():
    pass

interstellar = CustomInterstellar.from_SWG(args.attenuation_path,
                                    bounds = dict(column_density = (0.0,10.0)))

#absolute_alpha = xpsi.Parameter('absolute_alpha',
#                           strict_bounds = (0.5,1.5),
#                           bounds = (None, None),
#                           doc=('Absolute energy-independent scaling factor'),
#                           symbol = r'$\alpha$',
#                           value = None,
#                           permit_prepend = False)

signals = [[],]

alpha_bounds = dict(alpha = (0.1, 1.9))

#-------#
# NICER #
#-------#
if args.NICER:
    NICER = namespace()

    try:
        counts = np.loadtxt(args.NICER_matrix_path, dtype=np.double)
    except IOError:
        NICER.data = xpsi.Data.phase_bin__event_list(args.NICER_event_path,
                                              channels=np.arange(30, 150),
                                              phases=np.linspace(0.0, 1.0, 33),
                                              channel_column=1,
                                              phase_column=2,
                                              skiprows=3,
                                              dtype=np.double,
                                              first=0,
                                              last=119,
                                              exposure_time=1636760.0)

        np.savetxt(args.NICER_matrix_path, NICER.data.counts)
    else:
        NICER.data = xpsi.Data(counts,
                               channels=np.arange(30, 150),
                               phases=np.linspace(0.0, 1.0, 33),
                               first=0,
                               last=119,
                               exposure_time=1636760.0)

    NICER.instrument = CustomInstrument.NICER_XTI(bounds = alpha_bounds,
                                          values = {},
                                          ARF = args.NICER_arf_path,
                                          RMF = args.NICER_rmf_path,
                                          max_input = 1250,
                                          max_channel = 150,
                                          min_input = 0,
                                          min_channel = 30,
                                          channel_edges = args.NICER_channels_path,
                                          prefix = 'XTI')
    #NICER.instrument.merge(absolute_alpha)

    NICER.signal = CustomSignal(data = NICER.data,
                                  instrument = NICER.instrument,
                                  interstellar = interstellar,
                                  cache = False,
                                  workspace_intervals = 1000,
                                  epsrel = 1.0e-8,
                                  epsilon = 1.0e-3,
                                  sigmas = 10.0)

    signals[0].append(NICER.signal)

args.XMM = False
#-----#
# XMM #
#-----#
if args.XMM:
    def handle_XMM_event_list(event_list_path, instrument):
        """ Event list with channel given in eV in the last column. """

        events = np.loadtxt(event_list_path, dtype=np.int,
                              skiprows=3, usecols=-1)

        spectrum = np.zeros(len(instrument.channels), dtype=np.double)

        for event in events:
            for i in range(len(instrument.channel_edges) - 1):
                if instrument.channel_edges[i] <= event/1.0e3 < instrument.channel_edges[i+1]:
                    spectrum[i] += 1.0

        return spectrum.reshape(-1,1) # count spectrum

    #--------#
    # XMM-PN #
    #--------#
    PN = namespace()
    PN.instrument = CustomInstrument.XMM_PN(bounds = alpha_bounds,
                                            values = {},
                                            ARF = args.PN_arf_path,
                                            RMF = args.PN_rmf_path,
                                            max_input = 1194,
                                            max_channel = 299,
                                            min_input = 0,
                                            min_channel = 57,
                                            channel_edges = args.PN_channels_path,
                                            prefix = 'PN')
    #PN.instrument.merge(absolute_alpha)

    PN.data = xpsi.Data(handle_XMM_event_list(args.PN_spectrum_path, PN.instrument),
                        channels=PN.instrument.channels,
                        phases=np.array([0.0, 1.0]),
                        first=0,
                        last=len(PN.instrument.channels) - 1,
                        exposure_time=6.80873046875e3)

    spectrum = np.loadtxt(args.PN_background_path,
                          skiprows=3,
                          usecols=1,
                          dtype=np.double)[57:299]

    support = np.zeros((len(spectrum), 2), dtype=np.double)
    support[:,0] = spectrum - 4.0 * np.sqrt(spectrum)
    support[support[:,0] < 0.0, 0] = 0.0
    support[:,1] = spectrum + 4.0 * np.sqrt(spectrum)

    for i in range(support.shape[0]):
        if support[i,1] == 0.0:
            for j in range(i, support.shape[0]):
                if support[j,1] > 0.0:
                    support[i,0] = support[j,1]
                    break

    support *= 0.9212 * (PN.data.exposure_time / 4.51098e5) # BACKSCAL x exposure ratio

    support /= PN.data.exposure_time # need count rate, so divide by exposure time

    PN.signal = CustomSignal(data = PN.data,
                             instrument = PN.instrument,
                             interstellar = interstellar,
                             support = support,
                             cache = False,
                             workspace_intervals = 1000,
                             epsrel = 1.0e-8,
                             epsilon = 1.0e-3,
                             sigmas = 10.0)

    signals[0].append(PN.signal)

    #----------#
    # XMM-MOS1 #
    #----------#
    MOS1 = namespace()

    class derive(Derive):
        global PN

        def __init__(self):
            pass

        def __call__(self, boundto, caller=None):
            return PN.instrument['alpha']

    MOS1.instrument = CustomInstrument.XMM_MOS1(bounds = {'alpha': None},
                                            values = {'alpha': derive()},
                                            ARF = args.MOS1_arf_path,
                                            RMF = args.MOS1_rmf_path,
                                            max_input = 594,
                                            max_channel = 100,
                                            min_input = 2, # skip intervals 7 & 8
                                            min_channel = 20,
                                            channel_edges = args.MOS1_channels_path,
                                            prefix = 'MOS1')
    #MOS1.instrument.merge(absolute_alpha)

    MOS1.data = xpsi.Data(handle_XMM_event_list(args.MOS1_spectrum_path,
                                                MOS1.instrument),
                          channels=MOS1.instrument.channels,
                          phases=np.array([0.0, 1.0]),
                          first=0,
                          last=len(MOS1.instrument.channels) - 1,
                          exposure_time=1.795957421875e4)

    spectrum = np.loadtxt(args.MOS1_background_path,
                          skiprows=3,
                          usecols=1,
                          dtype=np.double)[20:100]

    support = np.zeros((len(spectrum), 2), dtype=np.double)
    support[:,0] = spectrum - 4.0 * np.sqrt(spectrum)
    support[support[:,0] < 0.0, 0] = 0.0
    support[:,1] = spectrum + 4.0 * np.sqrt(spectrum)

    for i in range(support.shape[0]):
        if support[i,1] == 0.0:
            for j in range(i, support.shape[0]):
                if support[j,1] > 0.0:
                    support[i,0] = support[j,1]
                    break

    support *= 1.074 * (MOS1.data.exposure_time / 1.57623e6) # BACKSCAL x exposure ratio

    support /= MOS1.data.exposure_time # need count rate, so divide by exposure time

    MOS1.signal = CustomSignal(data = MOS1.data,
                               instrument = MOS1.instrument,
                               interstellar = interstellar,
                               support = support,
                               cache = False,
                               workspace_intervals = 1000,
                               epsrel = 1.0e-8,
                               epsilon = 1.0e-3,
                               sigmas = 10.0)

    signals[0].append(MOS1.signal)

    #----------#
    # XMM-MOS2 #
    #----------#
    MOS2 = namespace()

    MOS2.instrument = CustomInstrument.XMM_MOS2(bounds = {'alpha': None},
                                            values = {'alpha': derive()},
                                            ARF = args.MOS2_arf_path,
                                            RMF = args.MOS2_rmf_path,
                                            max_input = 594,
                                            max_channel = 100,
                                            min_input = 2, # skip intervals 7 & 8
                                            min_channel = 20,
                                            channel_edges = args.MOS2_channels_path,
                                            prefix = 'MOS2')
    #MOS2.instrument.merge(absolute_alpha)

    MOS2.data = xpsi.Data(handle_XMM_event_list(args.MOS2_spectrum_path,
                                                MOS2.instrument),
                          channels=MOS2.instrument.channels,
                          phases=np.array([0.0, 1.0]),
                          first=0,
                          last=len(MOS2.instrument.channels) - 1,
                          exposure_time=1.8680734375e4)

    spectrum = np.loadtxt(args.MOS2_background_path,
                          skiprows=3,
                          usecols=1,
                          dtype=np.double)[20:100]

    support = np.zeros((len(spectrum), 2), dtype=np.double)
    support[:,0] = spectrum - 4.0 * np.sqrt(spectrum)
    support[support[:,0] < 0.0, 0] = 0.0
    support[:,1] = spectrum + 4.0 * np.sqrt(spectrum)

    for i in range(support.shape[0]):
        if support[i,1] == 0.0:
            for j in range(i, support.shape[0]):
                if support[j,1] > 0.0:
                    support[i,0] = support[j,1]
                    break

    support *= 1.260 * (MOS2.data.exposure_time / 1.51256e6) # BACKSCAL x exposure ratio

    support /= MOS2.data.exposure_time # need count rate, so divide by exposure time

    MOS2.signal = CustomSignal(data = MOS2.data,
                               instrument = MOS2.instrument,
                               interstellar = interstellar,
                               support = support,
                               cache = False,
                               workspace_intervals = 1000,
                               epsrel = 1.0e-8,
                               epsilon = 1.0e-3,
                               sigmas = 10.0)

    signals[0].append(MOS2.signal)

#-------#
# J0740 #
#-------#

bounds = dict(mass = (None, None),
              radius = (3.0*gravradius(1.0), 16.0),
              distance = (None, None),
              cos_inclination = (None, None))

spacetime = xpsi.Spacetime(bounds, dict(frequency = 346.53637))

bounds = dict(super_colatitude = (0.001, math.pi - 0.001),
              super_radius = (0.001, math.pi/2.0 - 0.001),
              phase_shift = (-0.25, 0.75),
              super_temperature = (5.1, 6.8))

primary = xpsi.HotRegion(bounds=bounds,
                            values={},
                            symmetry=True,
                            omit=False,
                            cede=False,
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=10,
                            max_sqrt_num_cells=64,
                            num_leaves=64,
                            num_rays=512,
                            is_secondary=False,
                            image_order_limit=3,
                            prefix='p')

bounds = dict(super_colatitude = (0.001, math.pi - 0.001),
              super_radius = (0.001, math.pi/2.0 - 0.001),
              phase_shift = (-0.25, 0.75),
              super_temperature = (5.1, 6.8))

secondary = xpsi.HotRegion(bounds=bounds,
                            values={},
                            symmetry=True,
                            omit=False,
                            cede=False,
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=10,
                            max_sqrt_num_cells=64,
                            num_leaves=64,
                            num_rays=512,
                            is_antiphased=True,
                            image_order_limit=3,
                            prefix='s')

from xpsi import HotRegions

hot = HotRegions((primary, secondary))

photosphere = CustomPhotosphere(hot = hot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))
photosphere.hot_atmosphere = args.atmosphere_path

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

prior = CustomPrior()

likelihood = xpsi.Likelihood(star = star, signals = signals,
                             num_energies = 128,
                             threads = 1,
                             externally_updated = True,
                             prior = prior)

# remember: get a point to check the likelihood between before using cluster
# or supercomp time
if args.XMM and args.NICER:
    p = [ 2.05908899, 11.68789646,  1.04426186,  0.04061792, -0.246576  ,
          0.66016391,  0.1450841 ,  6.02927387,  0.69387331,  1.54402738,
          0.10217615,  6.04285115,  1.12881329,  0.18888585,  1.0        ]
    likelihood.check(None, [-20698.938114135275], 1.0e-6,
                     physical_points=[p])
elif args.XMM:
    p = [ 2.05908899, 11.68789646,  1.04426186,  0.04061792, -0.246576  ,
          0.66016391,  0.1450841 ,  6.02927387,  0.69387331,  1.54402738,
          0.10217615,  6.04285115,  1.0,  0.18888585]
    likelihood.check(None, [-4563.4538676544], 1.0e-6,
                     physical_points=[p])
elif args.NICER:
    p = [ 2.05908899, 11.68789646,  1.04426186,  0.04061792, -0.246576  ,
          0.66016391,  0.1450841 ,  6.02927387,  0.69387331,  1.54402738,
          0.10217615,  6.04285115,  1.12881329,  0.18888585]
    likelihood.check(None, [-16135.484274586], 1.0e-6,
                     physical_points=[p])

if __name__ == '__main__': # sample from the posterior
    # transform relevant input information below to conmmand line arguments
    # and config file arguments

    wrapped_params = [0] * len(likelihood)
    wrapped_params[likelihood.index('p__phase_shift')] = 1
    wrapped_params[likelihood.index('s__phase_shift')] = 1

    if args.multinest:
        runtime_params = {'resume': args.resume,
                          'importance_nested_sampling': False,
                          'multimodal': False,
                          'n_clustering_params': None,
                          'outputfiles_basename': './run1_nlive1000_eff0.1_noCONST_noMM_noIS_tol-1',
                          'n_iter_before_update': 100,
                          'n_live_points': 1000,
                          'sampling_efficiency': 0.1,
                          'const_efficiency_mode': False,
                          'wrapped_params': wrapped_params,
                          'evidence_tolerance': 0.1,
                          'max_iter': -1,
                          'verbose': True}

        xpsi.Sample.nested(likelihood, prior, **runtime_params)

    elif args.emcee:
        try:
            backend
        except NameError:
            import emcee
            backend = emcee.backends.HDFBackend('samples.h5')

        chain = backend.get_chain()
        p = chain.reshape(-1,13)[np.argmax(backend.get_log_prob())]
        xpsi.ParameterSubspace.__call__(likelihood, p)
        print('Fiducial log-likelihood: %.6e' % likelihood())

        std = np.array(likelihood.vector) * 0.005

        runtime_params = {'resume': False,
                          'root_dir': './',
                          'nwalkers': 100,
                          'nsteps': 2000,
                          'walker_dist_moments': zip(likelihood.vector, std)}

        likelihood.threads = 1
        # Use MPI=False for testing purposes
        xpsi.MPI.COMM_WORLD.Barrier()
        backend = xpsi.Sample.ensemble(likelihood, prior,
                                       MPI=True, **runtime_params)
