//------------------------------------------------------------------------------
// Testing
//------------------------------------------------------------------------------

#include "house.h"

void house_run(HDIST dxi, HDIST dw, HDIST dn, HMOD* model, size_t K,
		double* xtrue, double* ztrue, const char* method_name,
		const char* file_name) {

	// Allocate for distributions
        HDIST* dxp = malloc(K * sizeof(HDIST));
        HDIST* dxu = malloc(K * sizeof(HDIST));

        // Timer variables
        clock_t t1, t2;

        // Run & time filter
        t1 = clock();
        for (size_t k = 0; k < K; k++) {
                if (k == 0)
                        dxp[k] = dxi;
                else
                        dxp[k] = predict(dxu[k-1], dw, k, model);
                dxu[k] = update(&ztrue[k*model->dim_measr], dxp[k],
                                dn, k, model);
        }
        t2 = clock();

	// Run time
	double runtime = (0.0 + t2 - t1) / CLOCKS_PER_SEC;

	// Print results to file
        print_dist(xtrue, dxu, model, K, runtime, method_name, file_name);

	// Free allocations
        for (size_t k = 1; k < K; k++) {
                free(dxp[k].mean);
                free(dxp[k].ccov);
                free(dxp[k].skew);
                free(dxp[k].kurt);
        }
        for (size_t k = 0; k < K; k++) {
                free(dxu[k].mean);
                free(dxu[k].ccov);
                free(dxu[k].skew);
                free(dxu[k].kurt);
        }
	free(dxu);
	free(dxp);

}

void house_test(HDIST dxi, HDIST dw, HDIST dn, HMOD* model, size_t K,
		double* xi, size_t num_methods, int* method,
		const char** method_names, const char** file_names) {

	// Generate true states and measurements
	double* xtrue = malloc(K * model->dim_state * sizeof(double));
	double* ztrue = malloc(K * model->dim_measr * sizeof(double));
	double* w;
	double* n;
	for (size_t k = 0; k < K; k++) {
		if (k == 0) {
			cblas_dcopy(model->dim_state, xi, 1, xtrue, 1);
		} else {
			w = rand_vect(dw);
			propagate(&xtrue[(k-1)*dxi.dim], &xtrue[k*dxi.dim],
				w, k, model);
			free(w);
		}
		n = rand_vect(dn);
		model->measurement(model->t0 + k * model->dt, &xtrue[k*dxi.dim],
			 n, model->params, &ztrue[k*model->dim_measr]);
		free(n);
	}

	// Run filter using each method
	for (size_t i = 0; i < num_methods; i++) {
		model->sig_method = method[i];
		house_run(dxi, dw, dn, model, K, xtrue, ztrue, method_names[i],
			file_names[i]);
	}

	// Free allocations
	free(xtrue);
	free(ztrue);

}

#define FORMAT "%15.6E "

// Save results of filter run to file
void print_dist(double* xtrue, HDIST* dx, HMOD* model, size_t K, double runtime,
		const char* method_name, const char* file_name) {

	// State dimension
	size_t N = model->dim_state;

	// Open file
	FILE* f = fopen(file_name, "w");

	// Header
	fprintf(f, "# Method: %s\n", method_name);
	fprintf(f, "# Run time: %12.6f s\n", runtime);
	int col = 1;
	fprintf(f, "# %3d: Time\n", col);
	col++;
	for (int i = 0; i < N; i++) {
		fprintf(f, "# %3d: True x%d\n", col, i+1);
		col++;
	}
	for (int i = 0; i < N; i++) {
                fprintf(f, "# %3d: Mean x%d\n", col, i+1);
                col++;
        }
	for (int i = 0; i < N; i++) {
                fprintf(f, "# %3d: StD  x%d\n", col, i+1);
                col++;
        }
	for (int i = 0; i < N; i++) {
                fprintf(f, "# %3d: Skew x%d\n", col, i+1);
                col++;
        }
	for (int i = 0; i < N; i++) {
                fprintf(f, "# %3d: Kurt x%d\n", col, i+1);
                col++;
        }
	/*
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
                	fprintf(f, "# %3d: Corr x%d x%d\n", col, i+1, j+1);
                	col++;
		}
        }
	*/

	// Print results
	double entry;
	for (int k = 0; k < K; k++) {

		// Times
		entry = model->t0 + k * model->dt;
		fprintf(f, FORMAT, entry);

		// True states
		for (size_t i = 0; i < N; i++) {
			entry = xtrue[N*k+i];
			fprintf(f, FORMAT, entry);
		}

		// Means
		for (size_t i = 0; i < N; i++) {
			entry = dx[k].mean[i];
			fprintf(f, FORMAT, entry);
		}

		// Standard deviations
		size_t ind = 0;
		for (size_t i = 0; i < N; i++) {
                        entry = dx[k].ccov[ind];
                        fprintf(f, FORMAT, entry);
			ind += N - i;
                }

		// Skewness
		for (size_t i = 0; i < N; i++) {
                        entry = dx[k].skew[i];
                        fprintf(f, FORMAT, entry);
                }

		// Kurtoses
		for (size_t i = 0; i < N; i++) {
                        entry = dx[k].kurt[i];
                        fprintf(f, FORMAT, entry);
                }

		/*
		// Correlations
		for (size_t i = 0; i < N; i++) {
			for (size_t j = 0; j < i; j++) {
				entry = covr[i][j] /
					sqrt(covr[i][i] * covr[j][j]);
				fprintf(f, FORMAT, entry);
			}
		}
		*/

		// Newline
		fprintf(f, "\n");

	}

	// Close file
	fclose(f);

}


