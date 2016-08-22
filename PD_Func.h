#ifndef PD_FUNC_H
#define PD_FUNC_H

// Compute the partial derivative of a vector valued function of position
// i.e. compute \partial f /  \partial x_{i} along the direction of x_{i}
// f is a mapping from R^{n} into R

namespace der_funcs{
	double partial_derivative(double (*the_func)(double *x, int n), double *x, double dx, int n, int dirn); // compute the value of the derivative of f at the point x by Richardson extrapolation
}

#endif 