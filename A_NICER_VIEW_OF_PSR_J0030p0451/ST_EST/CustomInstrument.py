from __future__ import print_function, division

import numpy as np
import math

import xpsi

class CustomInstrument(xpsi.Instrument):
    """ Methods and attributes specific to the NICER instrument.

    Currently tailored to the NICER light-curve SWG model specification.

    """
    def __init__(self, ratio, PI_channels, chan_edges, *args):
        """ Set channel edges attribute. """
        super(CustomInstrument, self).__init__(*args)
        self._ratio = ratio
        self._PI_channels = PI_channels
        self._chan_edges = chan_edges

        self._modified = self.matrix.copy()
        for i in range(self._modified.shape[0]):
            self._modified[i,:] *= self._ratio[i]

    @property
    def channels(self):
        return self._PI_channels

    @property
    def channel_edges(self):
        """ Get the channel edges. """
        return self._chan_edges

    def _construct_matrix(self, p):
        """ Implement response matrix parameterisation. """
        matrix = p[0]*p[1]*self._modified + (1.0 - p[1])*p[2]*self.matrix

        matrix[matrix < 0.0] = 0.0

        return matrix

    def __call__(self, p, signal, *args):
        """ Overwrite. """

        matrix = self._construct_matrix(p)

        self._folded_signal = np.dot(matrix, signal)

        return self._folded_signal

    @classmethod
    def from_SWG(cls, num_params, bounds,
                 ARF, RMF, ratio, max_input, min_input=0, chan_edges=None,
                 offset_correction=None):
        """ Constructor which converts files into :class:`numpy.ndarray`s.

        :param str ARF: Path to ARF which is compatible with
                                :func:`numpy.loadtxt`.

        :param str RMF: Path to RMF which is compatible with
                                :func:`numpy.loadtxt`.

        :param str ratio: Path to channel-by-channel ratio file.

        :param str chan_edges: Optional path to edges which is compatible with
                                :func:`numpy.loadtxt`.

        """
        try:
            ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
            RMF = np.loadtxt(RMF, dtype=np.double, skiprows=3, usecols=-1)
            ratio = np.loadtxt(ratio, dtype=np.double, skiprows=3)[:,2]
            if chan_edges:
                chan_edges = np.loadtxt(chan_edges, dtype=np.double, skiprows=3)
        except (OSError, IOError, TypeError, ValueError):
            print('A file could not be loaded.')
            raise

        matrix = np.zeros((1501,3980))

        for i in range(3980):
            matrix[:,i] = RMF[i*1501:(i+1)*1501]

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        edges = np.zeros(ARF[min_input:max_input,3].shape[0]+1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        RSP = np.ascontiguousarray(np.zeros(matrix[25:300,min_input:max_input].shape), dtype=np.double)

        for i in range(RSP.shape[0]):
            RSP[i,:] = matrix[i+25,min_input:max_input] * ARF[min_input:max_input,3] * 49.0/52.0

        PI_channels = np.arange(25, 300)

        ratios = ratio[:275]
        ratios[:10] = ratio[10]

        return cls(ratios, PI_channels, chan_edges[25:301,-2],
                   num_params, bounds, RSP, edges)
