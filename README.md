# SCOPE Shrinkage: A Unified Framework for Wavelet Denoising

## Overview
SCOPE (A Unified Framework for Wavelet Denoising) is a parametric wavelet shrinkage framework for signal denoising. It defines shrinkage rules using centered cumulative distribution functions (CDFs), yielding a flexible family of estimators:

$\delta (x;\lambda,k) = x |F^*(\lambda x)|^k$.

where $F^*(x) = 2F(x) - 1$, $\lambda > 0$ controls scale, and $k > 0$ controls transition sharpness.

By varying the underlying CDF (logistic, normal, Laplace, hyperbolic secant, Cauchy, uniform) and tuning the parameters ($\lambda$,$k$), SCOPE adapts to different signal structures and tail behaviors.

Parameters may be selected using oracle calibration or Stein’s unbiased risk estimate (SURE). Simulation studies show that SCOPE provides competitive performance across standard benchmark signals.

