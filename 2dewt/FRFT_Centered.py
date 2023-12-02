import numpy as np
import math
import os

def FRFT_Centered(x, alpha):

    '''
    inputs:
    x -> vector to be transformed
    alpha -> scaling factor

    outputs:
    y -> transformed vector
    '''
  
    x = np.array(x)
    N = len(x)

    factor_2 = np.exp(1j*np.pi*np.arange(N)*N*alpha)
    x_tilde = x*factor_2
    n = np.concatenate((np.arange(N), -np.arange(N, 0, -1)))
    factor = np.exp(-1j*np.pi*alpha*n**2)

    x_tilde = np.concatenate((x_tilde, np.zeros(N)))
    x_tilde = x_tilde * factor

    XX = np.fft.fft(x_tilde)
    YY = np.fft.fft(np.conj(factor))

    y = np.fft.ifft(XX * YY)
    y = y * factor

    y = y[:N]

    return y
