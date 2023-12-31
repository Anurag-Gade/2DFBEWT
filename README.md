# 2D fixed boundary point-based Empirical Wavelet Transform

**Code for the two-dimensional Empirical Wavelet Transform with fixed-boundary points (FBPs).**

There are various approaches to multiscale analysis, empirical wavelet transform being one. For multiscale analysis of signals, 1-dimension empirical wavelet transform (1DEWT) was introduced and for the analysis of images its two-dimensional counterpart was introduced. However in [1] and [2], two-dimensional empirical wavelet transform was used with **fixed boundary points (FBPs)**. Using frequency points say $f_p$, the FBPs are calculated by using the relation $B_p$ = 2 $\Pi$ $f_p$/N. 

## Instructions

The script to run is `main.py` which is included in the `src` folder. To run `main.py`, there are 4 arguments which are required;

- ***input_path***: Provide the absolute path to your input image.
- ***out_folder***: Provide the output folder path where the EWT modes will be saved.
- ***num_points***: Provide the number of boundary points.
- ***points***: Provide the frequency points $f_p$ in the form of a list i.e, [4,8,12,16,20]. 

## Citation
If this repository is useful to your research, please cite as below:

1. Gade, A., Dash, D. K., Kumari, T. M., Ghosh, S. K., Tripathy, R. K., & Pachori, R. B. (2023). Multiscale Analysis Domain Interpretable Deep Neural Network for Detection of Breast Cancer Using Thermogram Images. IEEE Transactions on Instrumentation and Measurement.

2. Muralidharan, N., Gupta, S., Prusty, M. R., & Tripathy, R. K. (2022). Detection of COVID19 from X-ray images using multiscale Deep Convolutional Neural Network. Applied Soft Computing, 119, 108610.

3. Gilles, J., Tran, G., & Osher, S. (2014). 2D empirical transforms. Wavelets, ridgelets, and curvelets revisited. SIAM Journal on Imaging Sciences, 7(1), 157-186.

4. Gilles, J. (2013). Empirical wavelet transform. IEEE transactions on signal processing, 61(16), 3999-4010.

