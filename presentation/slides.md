---
title: Quantifying the limits to population forecasts
subtitle: An example using the Yellowstone bison population
author: Andrew Tredennick
institute: Utah State University
date: \today
---

## Main goals

1. Fit Bayesian state-space model for bison population dynamics \alert{with} an environmental covariate
2. Compare out-of-sample forecasts with and without known environmental conditions
3. Partition forecast uncertainty into components

#	The data

##	Time series of bison counts (1970 - 2017)

Response
: Bison counts

Covariate
: Accumulated snow water equilivalent (West Yellowston SNOTEL)

\includegraphics[width=\textwidth]{../figures/bison_data_plots.png}

\begin{center}
\textcolor{blue}{training data} $\cdot$ \textcolor{orange}{validation data}
\end{center}

#	The model

##	Gompertz population growth

$$
\text{log}(z_{(t)}) \sim \text{Normal}\left( \text{log}(z_{(t-1)}) + r + b_0 \text{log}(z_{(t-1)}) + b_1 x_{(t)}, \sigma^2_\text{p} \right)
$$

$z_t$
: \alert{latent} population abundance in year *t*

$r$
: per capita growth rate

$b_0$
: density dependence

$b_1$
: effect of snow water equivalent

$x_t$
: accumulated snow water equivalent in year *t*

$\sigma^2_\text{p}$
: process variance

## Likelihood and fully specified model
 
 Likelihood

 $$
 y_{(t)} \sim \text{NB} \left(  z_{(t)} , \kappa \right)
 $$

 Full model

 $$
 \left[ \boldsymbol{\theta}_\text{p}, \kappa, z_{(t)}, z_{(t-1)} | y_{(t)}, x_{(t)} \right ] \propto \prod_{t=2}^{58} \underbrace{\left[ z_{(t)} | \boldsymbol{\theta}_\text{p}, z_{(t-1)}, x_{(t)} \right]}_{\text{process}} \prod_{t=1}^{48} \underbrace{\left[ y_{(t)} | \kappa, z_{(t)} \right]}_{\text{data}} \underbrace{\left[ \boldsymbol{\theta}_\text{p}, \kappa, z_{(t=1)}\right]}_{\text{parameters}}
 $$

 Includes a strong prior on $r$ based on Hobbs et al. 2015: $r \sim \text{Normal}(0.1, 0.02)$

# Results

## Posterior distributions of parameters

\includegraphics[width=\textwidth]{../figures/bison_post_params.png}

\alert{Note}: posterior distrbution of $r$ totally informed by prior.

## Posterior predictions, forecasts, and forecast partition

\includegraphics[height=\textheight]{../figures/bison_combined.png}

# Take home messages

## Take home

1. Snow water equivalent effect is weak -- right covariate?
2. Forecast uncertainty is large.
3. Forecast unceratinty dominated by (simulated) uncertainty of snow water equivalent


