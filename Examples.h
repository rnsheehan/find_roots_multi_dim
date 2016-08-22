#ifndef EXAMPLES_H
#define EXAMPLES_H

// Define the test functions for the problem

namespace examples{

	/* 2D Functions */

	double f11(double *x, int n);
	double f12(double *x, int n);

	double f21(double *x, int n);
	double f22(double *x, int n);

	double f31(double *x, int n);
	double f32(double *x, int n);

	double f41(double *x, int n);
	double f42(double *x, int n);

	/* 3D Functions */

	double g11(double *x, int n);
	double g12(double *x, int n);
	double g13(double *x, int n);

	double g21(double *x, int n);
	double g22(double *x, int n);
	double g23(double *x, int n);

	double g31(double *x, int n);
	double g32(double *x, int n);
	double g33(double *x, int n);

	double g41(double *x, int n);
	double g42(double *x, int n);
	double g43(double *x, int n);

	void test_partial_derivative(); 
	void test_coord_func_output(); 
	void test_jacobi_output(); 

	void test_2D1(); 
	void test_2D2(); 
	void test_2D3(); 
	void test_2D4(); 

	void test_3D1(); 
	void test_3D2(); 
	void test_3D3(); 
	void test_3D4(); 
}

#endif 