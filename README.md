# Shallow Water Model on a Yin–Yang Grid (McCormack Scheme)
# This repository contains the source code accompanying a manuscript currently under peer review.
## Overview

This repository contains a Fortran 90 implementation of a numerical model for solving the shallow water equations on a Yin–Yang grid. The model employs a fourth-order compact McCormack scheme combined with a stability enhancement strategy based on overlapping subdomains and inter-subdomain prognostic data transfer.

The implementation is designed for idealized global atmospheric dynamics experiments.

---

## Numerical Method

- Governing equations: Shallow Water Equations (SWE)
- Spherical geometry: Yin-Yang grid
- Spatial discretization: Fourth-order MacCormack (4/2) scheme
- Time integration: Fourth-order Runge Kutta scheme
- Stability strategy: Overlapping subdomain coupling with inter-subdomain prognostic variable exchange (Strategy 1)

---

## Test Cases

The model includes two standard benchmark test cases:

### 1. Rossby–Haurwitz Wave
A classical test for evaluating the model's ability to maintain large-scale wave structures on the sphere.

### 2. Midlatitude Unstable Tropospheric Jet 
An unstable tropospheric jet test case based on Gallewsky et al. (2004), used to evaluate the model's ability to capture dynamically unstable flow evolution.
## Requirements

- Fortran compiler (gfortran recommended)
- GNU Make
- Tecplot (for visualization of output files)

---

## Repository Structure

