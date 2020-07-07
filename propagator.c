#include "house.h"

// Propagate state
void propagate(double* xi, double* xf, double* w, size_t k, HMOD* model) {

	// Initial and final time
	double ti, tf;
	ti = model->t0 + k * model->dt;
	tf = ti + model->dt;

        // CVODE output
        int flag = 0;

        // Convert double arrays CVODE-friendly form
        N_Vector yi = N_VMake_Serial(model->dim_state, xi);
        N_Vector yf = N_VMake_Serial(model->dim_state, xf);
        N_Vector at = N_VMake_Serial(model->dim_state, model->abs_tol);

        // Choose linear multistep method
        int lmm = model->stiff ? CV_BDF : CV_ADAMS;

        // Nonlinear solver
        const int nls = CV_NEWTON;

        // Create CVODE object
        void* cvode_mem = CVodeCreate(lmm, nls);

        // Initialize CVODE solver
        flag = CVodeInit(cvode_mem, state_deriv_cvode, ti, yi);

        // Set noise & miscellaneous parameters
        HCPAR noisy_model;
        noisy_model.model = model;
        noisy_model.noise = w;
        flag = CVodeSetUserData(cvode_mem, &noisy_model);

        // Set tolerances
        flag = CVodeSVtolerances(cvode_mem, model->rel_tol, at);

	// Set maximum number of steps
	flag = CVodeSetMaxNumSteps(cvode_mem, model->max_steps);

        // Set up dense linear equation solver
        flag = CVDense(cvode_mem, model->dim_state);

        // Advance ODE solution
        double ts;
        flag = CVode(cvode_mem, tf, yf, &ts, CV_NORMAL);

	if (flag != 0) {
		printf("Error in CVODE solver\n");
		printf("x = %f, %f, %f, %f\n", xi[0], xi[1], xi[2], xi[3]);
		printf("w = %f\n", w[0]);
	}

        // Destroy CVODE object and vectors
        CVodeFree(&cvode_mem);
        N_VDestroy(at);
        N_VDestroy(yi);
        N_VDestroy(yf);

}

// CVODE-friendly state derivative function
int state_deriv_cvode(double t, N_Vector y, N_Vector yd, void* hcpar) {
	HCPAR hcp = *((HCPAR*) hcpar);
	double* x = N_VGetArrayPointer(y);
        double* v = N_VGetArrayPointer(yd);
        double* w = hcp.noise;
        void* p = hcp.model->params;
	return hcp.model->state_deriv(t, x, w, p, v);
}

