import numpy as np
import math
import os

def FRFT(X, alpha):

    X = X.flatten()
    N = len(X)

    n = np.concatenate((np.arange(0, N), np.arange(-N, 0)))
    factor = np.exp(-1j * np.pi * alpha * n**2)

    X_tilde = np.concatenate((X, np.zeros(N)))
    X_tilde = X_tilde * factor

    XX = np.fft.fft(X_tilde)
    YY = np.fft.fft(np.conj(factor))
    y = np.fft.ifft(XX * YY)
    y = y * factor
    y = y[:N]

    return y
