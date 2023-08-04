from __future__ import print_function, division

import numpy as np
import math

import xpsi

class CustomInterstellar(xpsi.Interstellar):
    """ Apply interstellar absorption. """

    def __init__(self, absorption, **kwargs):

        super(CustomInterstellar, self).__init__(**kwargs)

        self._supplied = absorption[0:351,:]

        self._energies = np.zeros(700, dtype=np.double)
        self._absorption = np.zeros(700, dtype=np.double)

        for i in range(self._supplied.shape[0]-1):
            att_diff = self._supplied[i+1, 1] - self._supplied[i, 1]
            E_diff = self._supplied[i+1, 0] - self._supplied[i, 0]
            self._absorption[2*i] = self._supplied[i,1] + 0.25*att_diff
            self._absorption[2*i+1] = self._supplied[i,1] + 0.75*att_diff
            self._energies[2*i] = self._supplied[i,0] + 0.25*E_diff
            self._energies[2*i+1] = self._supplied[i,0] + 0.75*E_diff

    @property
    def absorption(self):
        return self._absorption

    def __call__(self, p, channel_range, pulse):

        for i in range(pulse.shape[1]):
            pulse[:,i] *= self._absorption**(p[0]/0.4)

    @classmethod
    def from_SWG(cls, path, **kwargs):
        """ Load absorption file from the NICER SWG. """

        temp = np.loadtxt(path, dtype=np.double)

        absorption = temp[:,::2]

        return cls(absorption, **kwargs)
