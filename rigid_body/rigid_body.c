#include "house.h"

int f(double t, double* w, double* g, void* params, double* wd) {
        double* MOI = (double*) params;
        wd[0] = ((MOI[1] - MOI[2]) * w[1] * w[2] + g[0]) / MOI[0];
        wd[1] = ((MOI[2] - MOI[0]) * w[2] * w[0] + g[1]) / MOI[1];
        wd[2] = ((MOI[0] - MOI[1]) * w[0] * w[1] + g[2]) / MOI[2];
        return 0;
}

int h(double t, double* w, double* n, void* params, double* z) {
	z[0] = w[0] + n[0];
	return 0;
}

int main() {

	// Moments of inertia
        static double MOI[3];
        MOI[0] = 400;
        MOI[1] = 200;
        MOI[2] = 100;

	// Integration tolerances
        double reltol = 1E-12;
        static double abstol[3] = {1E-12, 1E-12, 1E-12};

        // System model
        HMOD M;
        M.dim_state = 3;          // State dimension
        M.dim_measr = 1;          // Measurement dimension
        M.dim_proc_noise = 3;     // Process noise dimension
        M.dim_meas_noise = 1;     // Measurement noise dimension
        M.state_deriv = &f;       // State derivative function
        M.measurement = &h;       // Measurement function
        M.t0 = 0;                 // Initial time
        M.dt = 0.1;               // Time step
        M.abs_tol = abstol;       // Absolute tolerance for ODE solver
        M.rel_tol = reltol;       // Relative tolerance for ODE solver
        M.stiff = 0;              // Indicates whether system is stiff
        M.sig_method = SIG_US;    // Method for generating sigma points
        M.params = MOI;           // Miscellaneous parameters
	M.max_steps = 1000; 	  // Maximum steps in ODE solver

	// Initial state distribution
	static double xim[3] = {0, 0, 0};
	static double xic[6] = {10, 0, 0, 10, 0, 10};
	static double xis[3] = {0.1, 0.1, 0.1};
	static double xik[3] = {4, 4, 4};
	HDIST dxi;
	dxi.dim = 3;
	dxi.mean = xim;
	dxi.ccov = xic;
	dxi.skew = xis;
	dxi.kurt = xik;

	// Measurement noise distribution
	static double nm = 0;
	static double nc = 0.01;
	static double ns = -0.1;
	static double nk = 200;
	HDIST dn;
	dn.dim = 1;
	dn.mean = &nm;
	dn.ccov = &nc;
	dn.skew = &ns;
	dn.kurt = &nk;

	// Process noise distribution
	static double wm[3] = {0, 0, 0};
        static double wc[6] = {1, 1, 1, 1, 1, 1};
        static double ws[3] = {0.1, 0.1, 0.1};
        static double wk[3] = {40, 40, 40};
        HDIST dw;
        dw.dim = 3;
        dw.mean = xim;
        dw.ccov = xic;
        dw.skew = xis;
        dw.kurt = xik;

	// Number of steps
        size_t K = 600;

        // Initial value of state
        static double xi[3] = {0, 2, 0};
        //double* xi = rand_vect(dxi);

        // Methods
        size_t num_methods = 2;
        static int method[2] = {SIG_HM, SIG_US};
        static const char* method_names[2] = {"HOUSE","UKF"};
        static const char* file_names[2] = {"out/dxhm", "out/dxus"};

        // Run test
        house_test(dxi, dw, dn, &M, K, xi, num_methods, method,
                method_names, file_names);

        //free(xi);

	return 0;

}
