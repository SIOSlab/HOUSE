#include "house.h"

// Parameter structure
typedef struct {
	double k;
	double l;
	double g;
	double mp;
	double mc;
} pend_const;

// State derivative function
int f(double t, double* x, double* w, void* params, double* v) {
	pend_const c = *((pend_const*) params);
	double xc, xcd, xcdd, th, thd, thdd, f_ext, sin_th, cos_th;
	xc = x[0];
	th = x[1];
	xcd = x[2];
	thd = x[3];
	f_ext = w[0];
	sin_th = sin(th);
	cos_th = cos(th);
	xcdd = (c.mp*(c.g*cos_th + c.l*thd*thd)*sin_th + f_ext - c.k*xc) /
		(c.mc + c.mp*sin_th*sin_th);
	thdd = - (xcdd*cos_th + c.g*sin_th) / c.l;
	/*
	xcdd = (-c.k*xc + f_ext-(c.mp*c.g*cos_th-c.mp*c.l*thd*thd)*sin_th) /
            (c.mc+c.mp*sin_th*sin_th);
        thdd = (c.g*sin_th - xcdd * cos_th) / c.l;
	*/
	v[0] = xcd;
	v[1] = thd;
	v[2] = xcdd;
	v[3] = thdd;
	return 0;
}

// Measurement function
int h(double t, double* x, double* n, void* params, double* z) {
	pend_const c = *((pend_const*) params);
	double xc, th, thd, xcdd, sin_th, cos_th;
        xc = x[0];
        th = x[1];
        thd = x[3];
        sin_th = sin(th);
        cos_th = cos(th);
	xcdd = (c.mp*(c.g*cos_th + c.l*thd*thd)*sin_th - c.k*xc) /
                (c.mc + c.mp*sin_th*sin_th);
	z[0] = xcdd;
	z[1] = sin_th;
	z[0] += n[0];
	z[1] += n[1];
	return 0;
}

int main() {

	// Pendulum constants
	static pend_const c;
	c.k = 1;
        c.l = 1;
        c.g = 9.81;
        c.mc = 2;
        c.mp = 0.5;

	// Integration tolerances
	double reltol = 1E-12;
	static double abstol[4] = {1E-12, 1E-12, 1E-12, 1E-12};

	// System model
	static HMOD M;
	M.dim_state = 4;          // State dimension
        M.dim_measr = 2;          // Measurement dimension
        M.dim_proc_noise = 1;     // Process noise dimension
        M.dim_meas_noise = 2;     // Measurement noise dimension
        M.state_deriv = &f;       // State derivative function
        M.measurement = &h;       // Measurement function
        M.t0 = 0;                 // Initial time
        M.dt = 0.1;               // Time step
        M.abs_tol = abstol;       // Absolute tolerance for ODE solver
        M.rel_tol = reltol;       // Relative tolerance for ODE solver
        M.stiff = 1;              // Indicates whether system is stiff
	M.max_steps = 10000;
        M.sig_method = SIG_US;    // Method for generating sigma points
        M.params = &c;	          // Miscellaneous parameters


	// Prior state distribution
	static double xim[4] = {0, 0, 0, 0};
	static double xic[10] = {10, 0, 0, 0,
				 10, 0, 0,
				 1, 0,
				 1};
	static double xis[4] = {0.1, 0.1, 0.1, 0.1};
	static double xik[4] = {4, 4, 4, 4};
	HDIST dxi;
	dxi.dim = 4;
	dxi.mean = xim;
	dxi.ccov = xic;
	dxi.skew = xis;
	dxi.kurt = xik;

	// Process noise distribution
	static double wm[1] = {0};
	static double wc[1] = {1};
	static double ws[1] = {-2};
	static double wk[1] = {20};
	HDIST dw;
	dw.dim = 1;
	dw.mean = wm;
	dw.ccov = wc;
	dw.skew = ws;
	dw.kurt = wk;

	// Measurement noise distribution
        static double nm[2] = {0, 0};
        static double nc[3] = {0.1, 0.1, 0.1};
        static double ns[2] = {1.5, -1};
        static double nk[2] = {10, 20};
        HDIST dn;
	dn.dim = 2;
        dn.mean = nm;
        dn.ccov = nc;
        dn.skew = ns;
        dn.kurt = nk;

	// Number of steps
	size_t K = 100;

	// Initial value of state
	//static double xi[4] = {0, 1, 0, 0};
	double* xi = rand_vect(dxi);

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
