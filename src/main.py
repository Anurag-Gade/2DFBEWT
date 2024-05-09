import numpy as np
import math 
import os
import cv2
import ast
import argparse
from EWT2D_LP_FilterBank import *
from PPFFT import *

parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description = "2D Empirical Wavelet Transform")

parser.add_argument('--input_path', action='store', required=True,
                    type=str, help='Path to the input file (image)')

parser.add_argument('--out_folder', action='store', required=True,
                    type=str, help='Path to store the output modes')

parser.add_argument('--num_points', action='store', required=True,
                    type=int, help='Number of boundary points')

parser.add_argument('--points', action='store', required=True,
                    type=str, help='Boundary points')

args = parser.parse_args()

input_path = args.input_path
out_folder = args.out_folder
num_points = args.num_points
points = args.points

points = ast.literal_eval(points)
points = [int(x) for x in points]

points = np.asarray(points)

f = cv2.imread(input_path)
f = cv2.cvtColor(f, cv2.COLOR_BGR2GRAY)

pseudoFFT = ppfft(f)
meanppfft = np.fft.fftshift(np.sum(np.abs(pseudoFFT), axis=1))
num = points * np.pi
points = num / np.round(len(meanppfft) / 2)

W = f.shape[1]
H = f.shape[0]

# Build the 2D filter bank
mfb = EWT2D_LP_FilterBank(points, H, W)

# We filter the signal to extract each subband
ff = np.fft.fft2(f)
ewtLP = [None] * len(mfb)
for k in range(len(mfb)-1):
    ewtLP[k] = np.real(np.fft.ifft2(np.conj(mfb[k]) * ff.T))

plt.figure()
for m in range(len(mfb)):
    # plt.subplot(2, 4, m+1)
    xx = ewtLP[m]
    # plt.imshow(xx)
    save_path = os.path.join(out_folder, "mode_{q}.png".format(q = m+1))
    cv2.imwrite(save_path, xx)

    print("Mode {num} saved!".format(num = m+1))








                    


