from __future__ import print_function, division

import numpy as np
import math

import xpsi

class CustomPhotosphere(xpsi.Photosphere):
    """ A photosphere extension to preload the numerical atmosphere NSX.

    """

    @xpsi.Photosphere.spot_atmosphere.setter
    def spot_atmosphere(self, path):
        NSX = np.loadtxt(path, dtype=np.double)
        logT = np.zeros(35)
        logg = np.zeros(11)
        mu = np.zeros(67)
        logE = np.zeros(166)

        reorder_buf = np.zeros((35,11,67,166))

        index = 0
        for i in range(reorder_buf.shape[0]):
            for j in range(reorder_buf.shape[1]):
                for k in range(reorder_buf.shape[3]):
                   for l in range(reorder_buf.shape[2]):
                        logT[i] = NSX[index,3]
                        logg[j] = NSX[index,4]
                        logE[k] = NSX[index,0]
                        mu[reorder_buf.shape[2] - l - 1] = NSX[index,1]
                        reorder_buf[i,j,reorder_buf.shape[2] - l - 1,k] = 10.0**(NSX[index,2])
                        index += 1

        buf = np.zeros(np.prod(reorder_buf.shape))

        bufdex = 0
        for i in range(reorder_buf.shape[0]):
            for j in range(reorder_buf.shape[1]):
                for k in range(reorder_buf.shape[2]):
                   for l in range(reorder_buf.shape[3]):
                        buf[bufdex] = reorder_buf[i,j,k,l]; bufdex += 1

        self._spot_atmosphere = (logT, logg, mu, logE, buf)
