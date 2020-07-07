#include "house.h"

#define G_0 9.80665

#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)

int f(double t, double* X, double* ad, void* params, double* Xd) {
	double B = *((double*) params);
	double* V = &X[3];
	Xd[0] = V[0];
	Xd[1] = V[1];
	Xd[2] = V[2];
	double v = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
	double* A = &Xd[3];
	A[0] =       - B * v * V[0] + ad[0];
	A[1] =       - B * v * V[1] + ad[1];
	A[2] = - G_0 - B * v * V[2] + ad[2];
        return 0;
}

#define THRESHOLD 1E-20

int h(double t, double* X, double* n, void* params, double* z) {
	double r = sqrt(X[0]*X[0] + X[1]*X[1]);
	if (r > THRESHOLD) {
		z[0] = atan2(X[1], -X[0]);
		z[1] = atan2(X[2], r);
	} else {
		z[0] = 0;
		z[1] = M_PI / 2;
	}
	z[0] += n[0];
	z[1] += n[1];
	return 0;
}

int main() {

	// Drag constant:  0.5 * A * CD * rho / m
	static double B = 1E-6;

	// Integration tolerances
        static double reltol = 1E-12;
        static double abstol[6] = {1E-12, 1E-12, 1E-12,
				   1E-12, 1E-12, 1E-12};

        // System model
        HMOD M;
        M.dim_state = 6;          // State dimension
        M.dim_measr = 2;          // Measurement dimension
        M.dim_proc_noise = 3;     // Process noise dimension
        M.dim_meas_noise = 2;     // Measurement noise dimension
        M.state_deriv = f;       // State derivative function
        M.measurement = h;       // Measurement function
        M.t0 = 0;                 // Initial time
        M.dt = 0.01; //0.01;               // Time step
        M.abs_tol = abstol;       // Absolute tolerance for ODE solver
        M.rel_tol = reltol;       // Relative tolerance for ODE solver
        M.stiff = 0;              // Indicates whether system is stiff
        M.sig_method = SIG_US;    // Method for generating sigma points
        M.params = &B;            // Miscellaneous parameters
	M.max_steps = 500; 	  // Maximum steps in ODE solver

	// Initial state distribution
	static double xim[6] = {0, 0, 0, 0, 0, 0};
	static double xic[21] = {1000, 0, 0, 0, 0, 0,
				 1000, 0, 0, 0, 0,
			   	 10, 0, 0, 0,
				 500,  0, 0,
				 500,  0,
				 1000};
	static double xis[6] = {0.1, 0.1, 0.1,
				0.1, 0.1, 0.1};
	static double xik[6] = {4, 4, 4,
				4, 4, 4,};
	HDIST dxi;
	dxi.dim = 6;
	dxi.mean = xim;
	dxi.ccov = xic;
	dxi.skew = xis;
	dxi.kurt = xik;

	// Measurement noise distribution
	static double nm[2] = {0, 0};
	static double nc[3] = {6 * ARC_MIN, 0, 6 * ARC_MIN};
	static double ns[2] = {-0.1, -0.1};
	static double nk[2] = {200, 200};
	HDIST dn;
	dn.dim = 2;
	dn.mean = nm;
	dn.ccov = nc;
	dn.skew = ns;
	dn.kurt = nk;

	// Process noise distribution
	static double wm[3] = {0, 0, 0};
        static double wc[6] = {0.1, 0, 0,
			       0.1, 0,
			       0.1};
        static double ws[3] = {0.1, 0.1, 0.1};
        static double wk[3] = {40, 40, 40};
        HDIST dw;
        dw.dim = 3;
        dw.mean = wm;
        dw.ccov = wc;
        dw.skew = ws;
        dw.kurt = wk;

	// Number of steps
        size_t K = 2500; // 2500;

        // Initial value of state
        double* xi = rand_vect(dxi);
	//static double xi[6] = {500, -500, 15, 100, 100, 200};

        // Methods
        size_t num_methods = 2;
        static int method[2] = {SIG_HM, SIG_US};
        static const char* method_names[2] = {"HOUSE","UKF"};
        static const char* file_names[2] = {"out/dxhm", "out/dxus"};

        // Run test
        house_test(dxi, dw, dn, &M, K, xi, num_methods, method,
                method_names, file_names);

        free(xi);

	return 0;

}
