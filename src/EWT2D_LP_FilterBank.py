import numpy as np
import math
import os
from EWT2D_LP_Scaling import *
from EWT2D_LP_Wavelet import *
from EWT2D_UP_LP_Wavelet import *


def EWT2D_LP_FilterBank(boundaries, W, H):
  
    Npic = len(boundaries)
  
    # Compute gamma
    gamma = 1
  
    for k in range(Npic-1):
        r = (boundaries[k+1] - boundaries[k]) / (boundaries[k+1] + boundaries[k])
        if r < gamma:
            gamma = r
          
    r = (math.pi - boundaries[Npic-1]) / (math.pi + boundaries[Npic-1])
    if r < gamma:
        gamma = r
      
    gamma = (1 - 1 / max(W, H)) * gamma  
  
    mfb = [None] * (Npic + 1)
    
    mfb[0] = EWT2D_LP_Scaling(boundaries[0], gamma, W, H)
  
    # Generate wavelets
  
    for k in range(Npic-1):
        mfb[k+1] = EWT2D_LP_Wavelet(boundaries[k], boundaries[k+1], gamma, W, H)
      
    mfb[Npic-1] = EWT2D_UP_LP_Wavelet(boundaries[Npic-1], gamma, W, H)
  
    return mfb

