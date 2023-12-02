# Littlewood-Paley Transform

import numpy as np
import math
import os

def EWT2D_LP_Wavelet(wn, wm, gamma, W, H):

    an = 1/(2*gamma*wn)
    am = 1/(2*gamma*wm)
    pbn = (1+gamma)*wn
    mbn = (1-gamma)*wn
    pbm = (1+gamma)*wm
    mbm = (1-gamma)*wm
    
    Mi = np.floor(W/2)
    Mj = np.floor(H/2)
    
    ymw = np.zeros((H,W))

    for i in range(0,W-1):
    
    for j in range(0,H-1):
    
      k1 = math.pi * (i-Mi)/Mi
      k2 = math.pi * (j-Mj)/Mj
      w = np.sqrt((k1**2) + (k2**2))
    
      if((w>=pbn) and (w<=mbm)):
    
          ymw[j+1, i+1] = 1
    
      elif((w>mbm) and (w<=pbm)):
    
          ymw[j+1, i+1] = np.cos(np.pi*EWT_beta(am*(w-mbm))/2)
    
      elif((w>=mbn) and (w<pbn)):
    
          ymw[j+1, i+1] = np.sin(np.pi*EWT_beta(an*(w-mbn))/2)
    
    ymw = np.fft.ifft(ymw)
    
    return ymw
