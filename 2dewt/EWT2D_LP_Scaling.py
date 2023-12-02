import numpy as np
import math
import os

def EWT2D_LP_Scaling(w1, gamma, W, H):
  
    an = 1/(2*gamma*w1)
    pbn = (1+gamma)*w1
    mbn = (1-gamma)*w1
  
    Mj = np.floor(H/2)
    Mi = np.floor(W/2)
    yms = np.zeros((H,W))
  
    for i in range(W):
      
        for j in range(H):
          
            k1 = np.pi*(i - Mi)/Mi
            k2 = np.pi*(j - Mj)/Mj
            w = np.sqrt(k1**2 + k2**2)
            if w < mbn:
                yms[j,i] = 1
            elif w >= mbn and w <= pbn:
                yms[j,i] = np.cos(np.pi * EWT_beta(an*(w-mbn)) / 2)
              
    yms = np.fft.ifftshift(yms)
  
    return yms

