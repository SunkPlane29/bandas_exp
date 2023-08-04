from __future__ import print_function, division

import numpy as np
import math

import xpsi

from xpsi import Parameter, make_verbose

class CustomInstrument(xpsi.Instrument):
    """ NICER XTI, and XMM pn, MOS1, and MOS2 instruments. """
    def construct_matrix(self):
        """ Implement response matrix parameterisation. """
        matrix = self['alpha'] * self.matrix # * self['absolute_alpha']
        matrix[matrix < 0.0] = 0.0

        return matrix

    def __call__(self, signal, *args):
        """ Overwrite. """

        matrix = self.construct_matrix()

        self._cached_signal = np.dot(matrix, signal)

        return self._cached_signal

    @classmethod
    @make_verbose('Loading NICER XTI response matrix',
                  'Response matrix loaded')
    def NICER_XTI(cls,
                  bounds, values,
                  ARF, RMF,
                  max_input,
                  max_channel,
                  min_input=0,
                  min_channel=0,
                  channel_edges=None,
                  **kwargs):
        """ Load NICER XTI instrument response matrix. """
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
        RMF = np.loadtxt(RMF, dtype=np.double, skiprows=3, usecols=-1)

        if channel_edges:
            channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

        matrix = np.zeros((1501,3451))

        for i in range(3451):
            matrix[:,i] = RMF[i*1501:(i+1)*1501]

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        RSP = np.zeros((max_channel - min_channel,
                        max_input - min_input), dtype=np.double)

        for i in range(RSP.shape[0]):
            RSP[i,:] = matrix[i+min_channel, min_input:max_input] * ARF[min_input:max_input,3] * 51.0/52.0

        channels = np.arange(min_channel, max_channel)

        alpha = Parameter('alpha',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('alpha', None),
                          doc='NICER XTI energy-independent scaling factor',
                          symbol = r'$\alpha_{\rm XTI}$',
                          value = values.get('alpha', None))

        return cls(RSP, edges, channels, channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)

    @classmethod
    @make_verbose('Loading XMM-pn response matrix',
                  'Response matrix loaded')
    def XMM_PN(cls,
               bounds, values,
               ARF, RMF,
               max_input,
               max_channel,
               min_input=0,
               min_channel=0,
               channel_edges=None,
               **kwargs):
        """ Load XMM-pn instrument response matrix. """
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
        RMF = np.loadtxt(RMF, dtype=np.double, skiprows=3, usecols=-1)

        if channel_edges:
            channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        matrix = np.zeros((4096, max_input))

        for i in range(max_input - min_input):
            matrix[:,i] = RMF[i*4096:(i+1)*4096]

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        RSP = np.zeros((max_channel - min_channel,
                        max_input - min_input), dtype=np.double)

        for i in range(RSP.shape[0]):
            RSP[i,:] = matrix[i+min_channel,min_input:max_input] * ARF[min_input:max_input,3]

        channels = np.arange(min_channel, max_channel)

        alpha = Parameter('alpha',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('alpha', None),
                          doc='XMM-pn energy-independent scaling factor',
                          symbol = r'$\alpha_{\rm pn}$',
                          value = values.get('alpha', None))

        return cls(RSP, edges, channels,
                   channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)

    @classmethod
    @make_verbose('Loading XMM-MOS1 response matrix',
                  'Response matrix loaded')
    def XMM_MOS1(cls,
                 bounds, values,
                 ARF, RMF,
                 max_input,
                 max_channel,
                 min_input=0,
                 min_channel=0,
                 channel_edges = None,
                 **kwargs):
        """ Load XMM-MOS1 instrument response matrix. """
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=9)
        RMF = np.loadtxt(RMF, dtype=np.double, skiprows=9, usecols=-1)

        if channel_edges:
            channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        matrix = np.zeros((800, max_input))

        for i in range(matrix.shape[1]):
            matrix[:,i] = RMF[i*800:(i+1)*800]

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        RSP = np.zeros((max_channel - min_channel,
                        max_input - min_input), dtype=np.double)

        for i in range(RSP.shape[0]):
            RSP[i,:] = matrix[i+min_channel,min_input:max_input] * ARF[min_input:max_input,3]

        channels = np.arange(min_channel, max_channel)

        alpha = Parameter('alpha',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('alpha', None),
                          doc='XMM-MOS1 energy-independent scaling factor',
                          symbol = r'$\alpha_{\rm MOS1}$',
                          value = values.get('alpha', None))

        return cls(RSP, edges, channels,
                   channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)

    @classmethod
    @make_verbose('Loading XMM-MOS2 response matrix',
                  'Response matrix loaded')
    def XMM_MOS2(cls,
                 bounds, values,
                 ARF, RMF,
                 max_input,
                 max_channel,
                 min_input=0,
                 min_channel=0,
                 channel_edges = None,
                 **kwargs):
        """ Load XMM-MOS2 instrument response matrix. """
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=9)
        RMF = np.loadtxt(RMF, dtype=np.double, skiprows=9, usecols=-1)

        if channel_edges:
            channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        matrix = np.zeros((800, max_input))

        for i in range(matrix.shape[1]):
            matrix[:,i] = RMF[i*800:(i+1)*800]

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        RSP = np.zeros((max_channel - min_channel,
                        max_input - min_input), dtype=np.double)

        for i in range(RSP.shape[0]):
            RSP[i,:] = matrix[i+min_channel,min_input:max_input] * ARF[min_input:max_input,3]

        channels = np.arange(min_channel, max_channel)

        alpha = Parameter('alpha',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('alpha', None),
                          doc='XMM-MOS2 energy-independent scaling factor',
                          symbol = r'$\alpha_{\rm MOS2}$',
                          value = values.get('alpha', None))

        return cls(RSP, edges, channels,
                   channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)
