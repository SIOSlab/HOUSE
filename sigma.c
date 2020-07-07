//------------------------------------------------------------------------------
// Sigma point generation
//------------------------------------------------------------------------------

#include "house.h"

// General sigma point generation
HSIG sigma_gen(HDIST state_dist, HDIST noise_dist, double* ax, double* bx,
		double* aw, double* bw) {
	size_t N, M, K;
        N = state_dist.dim;
        M = noise_dist.dim;
        K = 2 * (M + N) + 1;
        double* xs = calloc(K * N, sizeof(double));
        double* ws = calloc(K * M, sizeof(double));
	double* W = malloc(K * sizeof(double));
	double* col;
	double wtot = 0;
        for (size_t i = 0; i < N; i++) {
		col    = &xs[N * (i + 1)];
		col[i] =  ax[i];
		col    = &xs[N * (i + 1 + N + M)];
		col[i] = -bx[i];
		W[1 + i]         = 1 / (ax[i] * (ax[i] + bx[i]));
		W[1 + i + M + N] = 1 / (bx[i] * (ax[i] + bx[i]));
		wtot += 1 / (ax[i] * bx[i]);
	}
	for (size_t i = 0; i < M; i++) {
		col    = &ws[M * (i + 1 + N)];
		col[i] =  aw[i];
		col    = &ws[M * (i + K - M)];
		col[i] = -bw[i];
		W[1 + i + N]         = 1 / (aw[i] * (aw[i] + bw[i]));
		W[1 + i + N + M + N] = 1 / (bw[i] * (aw[i] + bw[i]));
		wtot += 1 / (aw[i] * bw[i]);
	}
	W[0] = 1 - wtot;
	HSIG sig;
	sig.dim_state = N;
	sig.dim_noise = M;
	sig.num_pts = K;
	sig.state = destandardize(xs, state_dist.mean, state_dist.ccov, K, N);
	sig.noise = destandardize(ws, noise_dist.mean, noise_dist.ccov, K, M);
	sig.weight = W;
	free(xs);
	free(ws);
	return sig;
}

// Generate sigma points for higher-moment unscented transform
HSIG sigma_hm(HDIST state_dist, HDIST noise_dist) {
	double k, s, ss, sk;

	double wtot = 0;
	for (size_t i = 0; i < state_dist.dim; i++) {
		s = state_dist.skew[i];
		k = state_dist.kurt[i];
		wtot += 1 / (k - s * s);
	}
	for (size_t i = 0; i < noise_dist.dim; i++) {
                s = noise_dist.skew[i];
                k = noise_dist.kurt[i];
                wtot += 1 / (k - s * s);
        }
	if (wtot > 1) {
		sk = wtot;
		ss = sqrt(sk);
	} else {
		sk = 1;
		ss = 1;
	}

	/*
	sk = 1;
	ss = 1;
	*/

	double* ax = malloc(state_dist.dim * sizeof(double));
	double* bx = malloc(state_dist.dim * sizeof(double));
	double* aw = malloc(noise_dist.dim * sizeof(double));
        double* bw = malloc(noise_dist.dim * sizeof(double));
	for (size_t i = 0; i < state_dist.dim; i++) {
		s = ss * state_dist.skew[i];
		k = sk * state_dist.kurt[i];
		ax[i] = (s + sqrt(4 * k - 3 * s * s)) / 2;
		bx[i] = ax[i] - s;
	}
	for (size_t i = 0; i < noise_dist.dim; i++) {
		s = ss * noise_dist.skew[i];
                k = sk * noise_dist.kurt[i];
                aw[i] = (s + sqrt(4 * k - 3 * s * s)) / 2;
                bw[i] = aw[i] - s;
	}
	HSIG sig = sigma_gen(state_dist, noise_dist, ax, bx, aw, bw);
	free(ax);
	free(bx);
	free(aw);
	free(bw);
	return sig;
}

// Generate sigma points for conventional unscented transform with n + kappa = 3
HSIG sigma_us(HDIST state_dist, HDIST noise_dist) {
        double* ax = malloc(state_dist.dim * sizeof(double));
        double* bx = malloc(state_dist.dim * sizeof(double));
        double* aw = malloc(noise_dist.dim * sizeof(double));
        double* bw = malloc(noise_dist.dim * sizeof(double));
	const double a = sqrt(3);
        for (size_t i = 0; i < state_dist.dim; i++) {
                ax[i] = a;
                bx[i] = a;
        }
        for (size_t i = 0; i < noise_dist.dim; i++) {
                aw[i] = a;
                bw[i] = a;
        }
        HSIG sig = sigma_gen(state_dist, noise_dist, ax, bx, aw, bw);
        free(ax);
        free(bx);
        free(aw);
        free(bw);
        return sig;
}

