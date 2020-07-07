//==============================================================================
// Matrix Statistics
//==============================================================================

#ifndef STATMAT_H
#define STATMAT_H

// Standard libraries
#include <stdio.h>
#include <stdlib.h>

// Linear algebra libraries
#include <cblas.h>
#include <lapacke.h>

// Mean
double* mean(double* x, double* w, size_t m, size_t n);

// Covariance (lower triangular)
double* cov(double* x, double* w, double* xm, size_t m, size_t n);

// Cholesky decomposition of covariance (lower triangular)
double* chol_cov(double* x, double* w, double* xm, size_t m, size_t n);

// Cross-covariance
double* cross_cov(double* x, double* y, double* w, double* xm, double* ym,
	size_t m, size_t nx, size_t ny);

// Linear minimum mean-square (LMMSE) estimate
double* lmmse(double* z, double* xm, double* zm, double* Pxz, double* Pzz,
	size_t m, size_t nx, size_t nz);

// Standardization
double* standardize(double* x, double* xm, double* Cxx, size_t m, size_t n);

// De-standardization
double* destandardize(double* x, double* xm, double* Cxx, size_t m, size_t n);

// Skewness
double* skewness(double* x, double* w, size_t m, size_t n);

// Kurtosis
double* kurtosis(double* x, double* w, size_t m, size_t n);

#endif
