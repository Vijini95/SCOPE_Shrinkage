# SCOPE Shrinkage: A Unified Framework for Wavelet Denoising

## Overview
SCOPE (A Unified Framework for Wavelet Denoising) is a parametric wavelet shrinkage framework for signal denoising. It defines shrinkage rules using centered cumulative distribution functions (CDFs), yielding a flexible family of estimators:

$\delta (x;\lambda,k) = x |F^*(\lambda x)|^k$,

where $F^*(x) = 2F(x) - 1$, $\lambda > 0$ controls scale, and $k > 0$ controls transition sharpness.

By varying the underlying CDF (logistic, normal, Laplace, hyperbolic secant, Cauchy, uniform) and tuning the parameters ($\lambda$, k), SCOPE adapts to different signal structures and tail behaviors.

Parameters may be selected using oracle calibration or Stein’s unbiased risk estimate (SURE). Simulation studies show that SCOPE provides competitive performance across standard benchmark signals.

## Methods
The SCOPE framework defines a parametric family of wavelet shrinkage rules based on centered cumulative distribution functions (CDFs). The method consists of the following components:

1. **Wavelet Transform:** The observed noisy signal is decomposed into multiple resolution levels using the discrete wavelet transform (DWT). This provides a multiscale representation in which noise and signal components can be separated in the coefficient domain.
   
2. **SCOPE Shrinkage Rule:** For each wavelet coefficient x, the SCOPE estimator is defined as
 
$\delta (x;\lambda,k) = x |F^*(\lambda x)|^k$,

where $F^*(x) = 2F(x) - 1$, $\lambda > 0$ controls scale, and $k > 0$ controls transition sharpness.

Different choices of F generate different shrinkage behaviors (logistic, normal, Laplace, hyperbolic secant, Cauchy, uniform).

3. **Parameter Calibration:** The shrinkage parameters ($\lambda$, k) can be selected using:
\item **Oracle calibration** (minimizing empirical MSE in simulations)
\item **Stein’s Unbiased Risk Estimate (SURE)** for data-driven tuning
The SURE criterion provides an unbiased estimate of the reconstruction risk under Gaussian noise.

4. **Signal Reconstruction:** After applying the shrinkage rule to the detail coefficients, the denoised signal is reconstructed using the inverse discrete wavelet transform (IDWT).
