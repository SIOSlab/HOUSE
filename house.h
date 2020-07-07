//==============================================================================
// Higher-Order UnScented Estimator (HOUSE)
//==============================================================================

#ifndef HOUSE_H
#define HOUSE_H

//------------------------------------------------------------------------------
// Libraries
//------------------------------------------------------------------------------

// Standard libraries
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// CBLAS & LAPACKE libraries for linear algebra
#include <cblas.h>
#include <lapacke.h>

// SUNDIALS & CVODE libraries for ODE solutions
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

// Matrix statistics
#include "statmat.h"

//------------------------------------------------------------------------------
// Structures & Type Definitions
//------------------------------------------------------------------------------

// Distribution
typedef struct {
	size_t dim;
	double* mean;
	double* ccov;
	double* skew;
	double* kurt;
} HDIST;

// Sigma points
typedef struct {
	size_t dim_state;
	size_t dim_noise;
	size_t num_pts;
	double* state;
	double* noise;
	double* weight;
} HSIG;

// Function of time, state, noise, and miscellaneous parameters
typedef int (*HFUN) (double t, double* state, double* noise, void* params,
	double* output);

// System model
typedef struct {
	int dim_state;		// State dimension
	int dim_measr;		// Measurement dimension
	int dim_proc_noise;	// Process noise dimension
	int dim_meas_noise;	// Measurement noise dimension
	HFUN state_deriv;	// State derivative function
	HFUN measurement; 	// Measurement function
	double t0;		// Initial time
	double dt;		// Time step
	double* abs_tol;	// Absolute tolerance for ODE solver
	double rel_tol;		// Relative tolerance for ODE solver
	int stiff;		// Indicates whether system is stiff
	int max_steps;		// Maximum number of steps in ODE solution
	int sig_method;		// Method for generating sigma points
	void* params;		// Miscellaneous parameters
} HMOD;

// Structure for passing process noise & system model to CVODE solver
typedef struct {
	HMOD* model;
	double* noise;
} HCPAR;

// Enumeration of sigma point generation methods
enum sig_method {SIG_US, SIG_HM};

//------------------------------------------------------------------------------
// Prediction & Update
//------------------------------------------------------------------------------

// State distribution prediction
HDIST predict(HDIST dxi, HDIST dw, size_t k, HMOD* model);

// State distribution update
HDIST update(double* Z, HDIST dxp, HDIST dn, size_t k, HMOD* model);

//------------------------------------------------------------------------------
// State Propagation
//------------------------------------------------------------------------------

// Propagate state
void propagate(double* xi, double* xf, double* w, size_t k, HMOD* model);

// CVODE-friendly state derivative function
int state_deriv_cvode(double t, N_Vector y, N_Vector yd, void* hcpar);

//------------------------------------------------------------------------------
// Sigma Point Generation
//------------------------------------------------------------------------------

// General sigma point generation
HSIG sigma_gen(HDIST state_dist, HDIST noise_dist, double* ax, double* bx,
                double* aw, double* bw);

// Generate sigma points for higher-moment unscented transform
HSIG sigma_hm(HDIST state_dist, HDIST noise_dist);

// Generate sigma points for conventional unscented transform with n + kappa = 3
HSIG sigma_us(HDIST state_dist, HDIST noise_dist);

//------------------------------------------------------------------------------
// Testing
//------------------------------------------------------------------------------

void house_run(HDIST dxi, HDIST dw, HDIST dn, HMOD* model, size_t K,
                double* xtrue, double* ztrue, const char* method_name,
                const char* file_name);

void house_test(HDIST dxi, HDIST dw, HDIST dn, HMOD* model, size_t K,
                double* xi, size_t num_methods, int* method,
                const char** method_names, const char** file_names);

void print_dist(double* xtrue, HDIST* dx, HMOD* model, size_t K, double runtime,
                const char* method_name, const char* file_name);

//-----------------------------------------------------------------------------
// Noisemaking
//-----------------------------------------------------------------------------

// Generates random vector from distribution
double* rand_vect(HDIST dist);

//------------------------------------------------------------------------------

#endif
