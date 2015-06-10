# Overview

This repository provides MATLAB code for estimating the parameters of the Generalized Quadratic Model (GQM), as described in the paper:

* IM Park, E Archer, N Priebe, JW Pillow, [Spectral methods for neural characterization using generalized quadratic models](http://pillowlab.princeton.edu/pubs/ParkI_GQM_NIPS2013.pdf), Advances in Neural Information Processing Systems 26, 2454--2462

The GQM is a statistical modeling framework arising from the neuroscience literature, and in particular the problem of modeling a neuron's response to a high-dimensional experimental stimulus. For instance: would you like to understand what a neuron in your visual cortex does when you watch random noise on your TV? If so, the GQM may be the right tool for you.

## The GQM

The GQM is a probabilistic model of a response $y$ conditioned upon a stimulus $x$, 

$$ y|x \sim P(f(Q(x)))$$

where $Q(x) = x^T C x + b^T x + a$ is a quadratic function of the stimulus and $P$ is a noise model. 

For example, for continuous-valued $y$ one might select $P$ to be Gaussian, whereas for discreted-valued $y$ one might choose $P$ to be Poisson. The GQM is closely-related to the Generalized Linear Model (GLM), the Linear-Nonlinear-Poisson model and the 2nd-order Volterra model, among others; please see the paper for further details. 


## Maximum Expected Log-Likelihood (MEL) Estimators 


The code in this repository implements several Maximum Expected Log Likelihood (MEL) estimators for the parameters $C$, $b$, and $a$. These estimators take different forms depending upon the distribution of $x$ and the noise model $P$.

MEL estimators provide a connection between the GQM and moment-based [dimensionality reduction](https://en.wikipedia.org/wiki/Dimensionality_reduction) techniques such as the [Spike-triggered covariance](https://en.wikipedia.org/wiki/Spike-triggered_covariance).



