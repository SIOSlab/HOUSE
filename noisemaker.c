#include "house.h"

//----------------------------------------------------------------------------
// Uniform & Normal Generators
//----------------------------------------------------------------------------

// Draws random number from uniform distribution on (0,1]
double uniform() {
        return (1.0 + rand()) / (1.0 + RAND_MAX);
}

// Draws random number from standard normal distribution
//      (uses Box-Muller transform)
double normal() {
        double u1, u2;
        u1 = uniform();
        u2 = uniform();
        return sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
}

//----------------------------------------------------------------------------
// Adapted from Joel Heinrich, "A Guide to the Pearson Type IV Distribution"
//	December 21, 2004 -- University of Pennsylvania
//	CDF/MEMO/STATISTICS/PUBLIC/6820
//----------------------------------------------------------------------------

double gammar2(double x,double y) {
	/* returns abs(gamma(x+iy)/gamma(x))^2 */
	const double y2=y*y, xmin = (2*y2>10.0) ? 2*y2 : 10.0;
	double r=1, s=1, p=1, f=0;
	while(x<xmin) {
		const double t = y/x++;
		r *= 1 + t*t;
	}
	while (p > s*DBL_EPSILON) {
		p *= y2 + f*f;
		p /= x++ * ++f;
		s += p;
	}
	return 1.0/(r*s);
}

double type4norm(double m,double nu,double a) {
	/* returns k */
	assert(m>0.5);
	return 0.5*M_2_SQRTPI*gammar2(m,0.5*nu)*exp(lgamma(m)-lgamma(m-0.5))/a;
}

double rpears4(double m,double nu,double a,double lam) {
	/* returns random Pearson IV deviate */
	const double k=type4norm(m,nu,1.0), b=2*m-2, M=atan(-nu/b);
	const double cosM=b/sqrt(b*b+nu*nu), r=b*log(cosM)-nu*M, rc=exp(-r)/k;
	double x,z;
	assert(m>1);
	do {
		int s=0;
		z = 0;
		if( (x=4*uniform()) > 2 ) {
			x -= 2;
			s = 1;
	}
	if (x > 1) x = 1 - (z=log(x-1)) ;
		x = (s) ? M + rc*x : M - rc*x;
	}
	while (fabs(x) >= M_PI_2 ||
		z + log(uniform()) > b*log(cos(x)) - nu*x - r);
	return a*tan(x) + lam;
}

//----------------------------------------------------------------------------
// Pearson Generator
//----------------------------------------------------------------------------

// Draws random number from Pearson distribution with mean 0 and variance 1
double pearson(const double skew, const double kurt) {
	double x, k, b1, b2;
	b1 = skew * skew;
	b2 = kurt;
	k = b1 * (b2 + 3) * (b2 + 3) /
		(4 * (4*b2 - 3*b1) * (2*b2 - 3*b1 - 6));
	if (k == 0) {
		if (b1 == 0 && b2 == 3) {
			x = normal();
		} else {
			printf("Distribution type not supported\n");
			exit(EXIT_FAILURE);
		}
	} else if (k > 0 && k < 1) {
		double r, m, nu, a, lam;
		r = 6 * (b2 - b1 - 1) / (2*b2 - 3*b1 - 6);
		m = r/2 + 1;
		nu = -r * (r - 2) * sqrt(b1) /
			sqrt(16*(r - 1) - b1*(r - 2)*(r - 2));
		a = sqrt(16*(r - 1) - b1*(r - 2)*(r - 2)) / 4;
		lam = -(r - 2) * sqrt(b1) / 4;
		x = rpears4(m, nu, a, lam);
	} else {
		printf("Distribution type not supported\n");
                exit(EXIT_FAILURE);
	}
	return x;
}

//-----------------------------------------------------------------------------
// Noisemaking
//-----------------------------------------------------------------------------

// Generates random vector from distribution
double* rand_vect(HDIST dist) {
	double* xs = malloc(dist.dim * sizeof(double));
	for (size_t i = 0; i < dist.dim; i++) {
		xs[i] = pearson(dist.skew[i], dist.kurt[i]);
	}
	double* x = destandardize(xs, dist.mean, dist.ccov, 1, dist.dim);
	free(xs);
	return x;
}

//-----------------------------------------------------------------------------
