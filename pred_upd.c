//------------------------------------------------------------------------------
// Prediction & Update
//------------------------------------------------------------------------------

#include "house.h"

// State distribution prediction
HDIST predict(HDIST dxi, HDIST dw, size_t k, HMOD* model) {

	// State & noise dimension
	size_t N = dxi.dim;
	size_t M = dw.dim;

	// Generate sigma points
	HSIG sig;
	switch (model->sig_method) {
		case SIG_HM:
			sig = sigma_hm(dxi, dw);
			break;
		case SIG_US:
			sig = sigma_us(dxi, dw);
			break;
		default:
			sig = sigma_us(dxi, dw);
	}

	// Compute predicted states at sigma points
	double* xp = malloc(N * sig.num_pts * sizeof(double));
	for (size_t i = 0; i < sig.num_pts; i++)
		propagate(&sig.state[N*i], &xp[N*i], &sig.noise[M*i], k, model);

	// Predicted state distribution
	HDIST dxp;
	dxp.dim = dxi.dim;

	// Predicted state mean
	dxp.mean = mean(xp, sig.weight, sig.num_pts, N);

	// Square root of predicted state covariance
	dxp.ccov = chol_cov(xp, sig.weight, dxp.mean, sig.num_pts, N);

	// Standardized predicted state
	double* xs = standardize(xp, dxp.mean, dxp.ccov, sig.num_pts, N);

	// Skewness and kurtosis of standardized predicted state
	dxp.skew = skewness(xs, sig.weight, sig.num_pts, N);
	dxp.kurt = kurtosis(xs, sig.weight, sig.num_pts, N);

	// Free allocations
	free(xp);
	free(xs);
	free(sig.state);
	free(sig.noise);
	free(sig.weight);

	// Return predicted state distribution
	return dxp;

}

// State distribution update
HDIST update(double* Z, HDIST dxp, HDIST dn, size_t k, HMOD* model) {

	// Time
	double t = model->t0 + k * model->dt;

	// State & noise dimension
        size_t N = dxp.dim;
        size_t M = dn.dim;
	size_t K = model->dim_measr;

        // Generate sigma points
        HSIG sig;
        switch (model->sig_method) {
                case SIG_HM:
                        sig = sigma_hm(dxp, dn);
                        break;
                case SIG_US:
                        sig = sigma_us(dxp, dn);
                        break;
                default:
                        sig = sigma_us(dxp, dn);
        }

	// Measurements at sigma points
	double* z = malloc(K * sig.num_pts * sizeof(double));
	int flag;
	for (size_t i = 0; i < sig.num_pts; i++) {
		flag = model->measurement(t, &sig.state[N*i], &sig.noise[M*i],
			model->params, &z[K*i]);
		if (flag != 0) printf("Error in measurement function\n");
	}

	// Measurement mean
	double* zm = mean(z, sig.weight, sig.num_pts, K);

	// Measurement covariance
	double* Pzz = cov(z, sig.weight, zm, sig.num_pts, K);

	// State-measurement cross-covariance
	double* Pxz = cross_cov(sig.state, z, sig.weight, dxp.mean, zm,
			sig.num_pts, N, K);

	// LMMSE states at sigma points
	double* xl = lmmse(z, dxp.mean, zm, Pxz, Pzz, sig.num_pts, N, K);

	// LMMSE error
	double* xe = malloc(N * sig.num_pts * sizeof(double));
	cblas_dcopy(N * sig.num_pts, sig.state, 1, xe, 1);
	cblas_daxpy(N * sig.num_pts, -1, xl, 1, xe, 1);

	// LMMSE error mean (should be zero)
	double* xem = mean(xe, sig.weight, sig.num_pts, N);

	// Updated state distribution
	HDIST dxu;
	dxu.dim = dxp.dim;

	// Square root of LMMSE state error covariance
	dxu.ccov = chol_cov(xe, sig.weight, xem, sig.num_pts, N);

	// Standardized LMMSE state error
	double* xes = standardize(xe, xem, dxu.ccov, sig.num_pts, N);

	// Skewness and covariance of standardized LMMSE state error
	dxu.skew = skewness(xes, sig.weight, sig.num_pts, N);
	dxu.kurt = kurtosis(xes, sig.weight, sig.num_pts, N);

	// Updated state mean
	dxu.mean = lmmse(Z, dxp.mean, zm, Pxz, Pzz, 1, N, K);

	// Free allocations
	free(z);
	free(zm);
	free(Pzz);
	free(Pxz);
	free(xl);
	free(xe);
	free(xem);
	free(xes);
	free(sig.state);
	free(sig.noise);
	free(sig.weight);

	// Return updated measurement distribution
	return dxu;

}

