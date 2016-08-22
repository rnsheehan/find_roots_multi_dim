#ifndef ROOT_MD_H
#define ROOT_MD_H

// Declaration of the class used to find the roots of multidimensional functions
// Given a function F(X), with F = { f_{1}(x, y, z), f_{2}(x, y, z), f_{3}(x, y, z) }^{T}, X = {x, y, z}
// Find the roots of the non-linear system of equations F(X) = 0
// Methods include Newton-Raphson, Broyden with Sherman-Morrison Inverse, Gradient Search and combination Gradient + Broyden Search
// R. Sheehan 25 - 3 - 2013

typedef double (*coordfunc)(double *x, int n); // this maps a vector x into a real number 

class find_root_multi_D{

public:
	// Constructor
	find_root_multi_D();
	find_root_multi_D(double (*func1)(double *, int ), double (*func2)(double *, int ));
	find_root_multi_D(double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int ));

	// Deconstructor
	//~find_root_multi_D(){ delete[] F; delete[] dx; delete[] xold; delete[] xnew; delete[] J; }

	// Methods
	void newton_raphson_search(double *xinit, double toler, bool loud = false); 
	void broyden_search(double *xinit, double toler, bool loud = false);
	void gradient_search(double *xinit, double toler, bool loud = false);
	void gradient_broyden_search(double *xinit, double toler, bool loud = false);

	// these functions can be declared private once you know they work correctly
	// end-user does not need access to these functions
//private:
	double sum_sqr(double *x); // evaluate the sum over the squares of the coordinate functions at x

	void fill_F(double *x, bool loud = false); // return the value of the function at the point x
	void jacobi_mat(double *pos, bool loud = false); // what is the jacobi matrix associated with the mapping
	void sherman_morrison_inverse(double **Aold, double **Anew, double *dxold, double *dF, bool loud = false); // estimate the inverse of a matrix using the Sherman-Morrison formula

private:
	int ndim; // size of search space 2 or 3 dimensions3
	int maxit; 

	double *F; 
	double *dx;
	double *xold; 
	double *xnew;
	double **J; 

	// the coordinate functions are stored as elements of the vector func_vec
	coordfunc *func_vec; 
};

#endif