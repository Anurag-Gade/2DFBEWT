import numpy as np
import math
import os
from EWT_beta import *

def EWT2D_UP_LP_Wavelet(wn, gamma, W, H):

    an=1/(2*gamma*wn)
    pbn=(1+gamma)*wn
    mbn=(1-gamma)*wn

    Mi = np.floor(W/2)
    Mj = np.floor(H/2)

    ymw = np.ones((H,W))

    for i in range(0, W-1):

        for j in range(0,H-1):

            k1 = np.pi*(i-Mi)/Mi
            k2 = np.pi*(j-Mj)/Mj

            w = np.sqrt((k1**2) + (k2**2))

            if(w<mbn):
                ymw[j+1,i+1] = 0

            elif((w>=mbn) and (w<pbn)):
                ymw[j+1,i+1] = np.sin(np.pi*EWT_beta(an*(w-mbn))/2)

    ymw = np.fft.ifftshift(ymw)

    return ymw
