import numpy as np
import math
import os
from FRFT import *
from FRFT_Centered import *

def ppfft(X,S1=1,S2=1):

    '''
    Initializing input and output array sizes
    '''
  
    input_shape = X.shape
    N1 = input_shape[0]
    N2 = input_shape[1]

    N = math.ceil(max(input_shape[0], input_shape[1])/2)*2
    Xnew = np.zeros((N,N))

    Xnew[int(N//2 - np.floor(N1/2)) : int(N//2 - np.floor(N1/2) + N1), int(N//2 - np.floor(N2/2)) : int(N//2 - np.floor(N2/2) + N2)] = X
    X = Xnew
    Y = np.zeros((2*S1*N, 2*S2*N))

    '''
    Constructing quadrant 1 and quadrant 3
    '''
  
    f_tilde = np.fft.fft(np.concatenate((X, np.zeros(((S1*2-1)*N, N)))), axis=0)
    f_tilde = np.fft.fftshift(f_tilde, axes=1)
    f_tilde = np.concatenate((f_tilde, np.zeros((2*S1*N, N*(S2-1)))), axis=1)

    for i in range(-N*S1, N*S1):
      
        Y[i+N*S1+1, S2*N-1::-1] = np.transpose(FRFT_Centered(f_tilde[i+N*S1+1, :], i/(N**2*S1*S2)))

    '''
    Constructing quadrant 2 and 4
    '''
    f_tilde = np.fft.fft(np.concatenate([X, np.zeros((N, (S1*2-1)*N))], axis=1), axis=1)
    f_tilde = np.fft.fftshift(f_tilde, axes=1)
    f_tilde = np.transpose(f_tilde)
    f_tilde = np.concatenate([f_tilde, np.zeros((2*S1*N, N*(S2-1)))], axis=1)

    Y = np.zeros((2*N*S1, 2*N*S2), dtype=np.complex128)
  
    for ll in range(-N*S1, N*S1):
      
        factor = np.exp(1j*2*np.pi*np.arange(S2*N)*(N*S2/2-1)*ll/(N**2*S1*S2))
        Y[ll+N*S1, S2*N:2*N*S2] = np.transpose(FRFT(f_tilde[ll+N*S1, :]*factor, ll/(N**2*S1*S2)))

    return Y
