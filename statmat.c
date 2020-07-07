//==============================================================================
// Matrix Statistics
//==============================================================================

#include "statmat.h"

// Mean
double* mean(double* x, double* w, size_t m, size_t n) {
	double* xm = calloc(n, sizeof(double));
	cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, 1, x, n, w, 1, 0, xm, 1);
	return xm;
}

// Covariance (lower triangular)
double* cov(double* x, double* w, double* xm, size_t m, size_t n) {
	double* Pxx = calloc(n * (n + 1) / 2, sizeof(double));
	for (size_t i = 0; i < m; i++)
		cblas_dspr(CblasColMajor, CblasLower, n,
			w[i], &x[i*n], 1, Pxx);
	cblas_dspr(CblasColMajor, CblasLower, n, -1, xm, 1, Pxx);
	return Pxx;
}

// Cholesky decomposition of covariance (lower triangular)
double* chol_cov(double* x, double* w, double* xm, size_t m, size_t n) {
	double* Cxx = cov(x, w, xm, m, n);
	int flag = LAPACKE_dpptrf(LAPACK_COL_MAJOR, 'L', n, Cxx);
	if (flag != 0) printf("Error in LAPACKE_dpptrf: %d\n", flag);
	return Cxx;
}

// Cross-covariance
double* cross_cov(double* x, double* y, double* w, double* xm, double* ym,
		size_t m, size_t nx, size_t ny) {
	double* Pxy = calloc(nx * ny, sizeof(double));
	for (size_t i = 0; i < m; i++)
		cblas_dger(CblasColMajor, nx, ny,
			w[i], &x[i*nx], 1, &y[i*ny], 1, Pxy, nx);
	cblas_dger(CblasColMajor, nx, ny, -1, xm, 1, ym, 1, Pxy, nx);
	return Pxy;
}

// Linear minimum mean-square (LMMSE) estimate
double* lmmse(double* z, double* xm, double* zm, double* Pxz, double* Pzz,
		size_t m, size_t nx, size_t nz) {
	double* xe = malloc(m * nx * sizeof(double));
	double* ze = malloc(nz * sizeof(double));
	double* Czz = malloc((nz * (nz + 1) / 2) * sizeof(double));
	int flag;
	for (size_t i = 0; i < m; i++) {
		cblas_dcopy(nz * (nz + 1) / 2, Pzz, 1, Czz, 1);
		cblas_dcopy(nx, xm, 1, &xe[i * nx], 1);
		cblas_dcopy(nz, &z[i * nz], 1, ze, 1);
		cblas_daxpy(nz, -1, zm, 1, ze, 1);
		flag = LAPACKE_dppsv(LAPACK_COL_MAJOR, 'L', nz, 1, Czz,
				ze, nz);
		if (flag != 0) printf("Error in LAPACKE_dppsv: %d\n", flag);
		cblas_dgemv(CblasColMajor, CblasNoTrans,
			nx, nz, 1, Pxz, nx, ze, 1, 1, &xe[i * nx], 1);
	}
	free(ze);
	free(Czz);
	return xe;
}

// Standardization
double* standardize(double* x, double* xm, double* Cxx, size_t m, size_t n) {
	double* xs = malloc(n * m * sizeof(double));
	cblas_dcopy(n * m, x, 1, xs, 1);
	for (size_t i = 0; i < m; i++) {
		cblas_daxpy(n, -1, xm, 1, &xs[i * n], 1);
		cblas_dtpsv(CblasColMajor, CblasLower, CblasNoTrans,
			CblasNonUnit, n, Cxx, &xs[i * n], 1);
	}
	return xs;
}

// De-standardization
double* destandardize(double* x, double* xm, double* Cxx, size_t m, size_t n) {
	double* xd = malloc(n * m * sizeof(double));
	cblas_dcopy(n * m, x, 1, xd, 1);
	for (size_t i = 0; i < m; i++) {
		cblas_dtpmv(CblasColMajor, CblasLower, CblasNoTrans,
			CblasNonUnit, n, Cxx, &xd[i * n], 1);
		cblas_daxpy(n, 1, xm, 1, &xd[i * n], 1);
	}
	return xd;
}

// Skewness
double* skewness(double* x, double* w, size_t m, size_t n) {
	double* x3 = malloc(m * n * sizeof(double));
	double* s = calloc(n, sizeof(double));
	cblas_dcopy(m * n, x, 1, x3, 1);
	cblas_dtbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		m * n, 0, x, 1, x3, 1);
	cblas_dtbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                m * n, 0, x, 1, x3, 1);
	cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, 1, x3, n, w, 1, 0, s, 1);
	free(x3);
	return s;
}

// Kurtosis
double* kurtosis(double* x, double* w, size_t m, size_t n) {
	double* x4 = malloc(m * n * sizeof(double));
        double* s = calloc(n, sizeof(double));
        cblas_dcopy(m * n, x, 1, x4, 1);
        cblas_dtbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                m * n, 0, x, 1, x4, 1);
        cblas_dtbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                m * n, 0, x, 1, x4, 1);
	cblas_dtbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                m * n, 0, x, 1, x4, 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, 1, x4, n, w, 1, 0, s, 1);
        free(x4);
        return s;
}
