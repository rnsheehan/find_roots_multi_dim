#ifndef ATTACH_H
#include "Attach.h"
#endif

/* 2D Function 1 */

double examples::f11(double *x, int n)
{
	// Coordinate function for the vector valued function F = (f11, f12)^{T}

	return 	( 3.0*template_funcs::DSQR(x[1]) - template_funcs::DSQR(x[2])); 
}

double examples::f12(double *x, int n)
{
	// Coordinate function for the vector valued function F = (f11, f12)^{T}

	return ( 3.0* x[1] * template_funcs::DSQR(x[2]) - x[1] * template_funcs::DSQR(x[1]) - 1.0 ); 
}

/* 2D Function 2 */

double examples::f21(double *x, int n)
{
	// Coordinate function for the vector valued function F = (f11, f12)^{T}

	return 	( log( template_funcs::DSQR(x[1]) + template_funcs::DSQR(x[2])) - sin(x[1]*x[2]) - log(2.0) - log(PI) ); 
}

double examples::f22(double *x, int n)
{
	// Coordinate function for the vector valued function F = (f11, f12)^{T}

	return ( exp( x[1] - x[2] ) + cos( x[1] * x[2] ) ); 
}

/* 2D Function 3 */

double examples::f31(double *x, int n)
{
	// Coordinate function for the vector valued function F = (f11, f12)^{T}

	return 	( x[1] * ( 1.0 - x[1] ) + 4.0 * x[2] - 12.0 ); 
}

double examples::f32(double *x, int n)
{
	// Coordinate function for the vector valued function F = (f11, f12)^{T}

	return ( template_funcs::DSQR( x[1] - 2.0 ) + template_funcs::DSQR( 2.0 * x[2] - 3.0 ) - 25.0); 
}

/* 2D Function 4 */

double examples::f41(double *x, int n)
{
	// Coordinate function for the vector valued function F = (f11, f12)^{T}

	return 	( 5.0*template_funcs::DSQR(x[1])-template_funcs::DSQR(x[2])); 
}

double examples::f42(double *x, int n)
{
	// Coordinate function for the vector valued function F = (f11, f12)^{T}

	return (x[2]-0.25*(sin(x[1])+cos(x[2]))); 
}

/* 3D Function 1 */

double examples::g11(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g11, g12, g13)^{T}

	return (template_funcs::DSQR(x[1])+x[2]-37.0);
}

double examples::g12(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g11, g12, g13)^{T}

	return (x[1]-template_funcs::DSQR(x[2])-5.0);
}

double examples::g13(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g11, g12, g13)^{T}

	return (x[1]+x[2]+x[3]-3.0);
}

/* 3D Function 2 */

double examples::g21(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g21, g22, g23)^{T}

	return (6.0*x[1]-2.0*cos(x[2]*x[3])-1.0);
}

double examples::g22(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g21, g22, g23)^{T}

	return (9.0*x[2]+sqrt(template_funcs::DSQR(x[1])+sin(x[3])+1.06)+0.9);
}

double examples::g23(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g21, g22, g23)^{T}

	return (60.0*x[3]+3.0*exp(-1.0*x[1]*x[2])-3.0+10.0*PI);
}

/* 3D Function 3 */

double examples::g31(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g31, g32, g33)^{T}

	return (15.0*x[1]+template_funcs::DSQR(x[2])-4.0*x[3]-13.0);
}

double examples::g32(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g31, g32, g33)^{T}

	return (template_funcs::DSQR(x[1])+10.0*x[2]-x[3]-11.0);
}

double examples::g33(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g31, g32, g33)^{T}

	return (x[2]*template_funcs::DSQR(x[2])-25.0*x[3]+22.0);
}

/* 3D Function 4 */

double examples::g41(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g41, g42, g43)^{T}

	return (10.0*x[1]-2.0*template_funcs::DSQR(x[2])+x[2]-2.0*x[3]-5.0);
}

double examples::g42(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g41, g42, g43)^{T}

	return (8.0*template_funcs::DSQR(x[2])+4.0*template_funcs::DSQR(x[3])-9.0);
}

double examples::g43(double *x, int n)
{
	// Coordinate function for the vector valued function F = (g41, g42, g43)^{T}

	return (8.0*x[2]*x[3]+4.0);
}

void examples::test_partial_derivative()
{
	// run a test on the partial derivative calculator to ensure it is running correctly

	double *x = new (double [2+1]); 

	x[1] = 7.0; x[2] = 13.0; 

	// \partial f_{11} / \partial x taken at x = 7, y = 13
	std::cout<< der_funcs::partial_derivative(*f11,x,0.1,2,1) <<"\n"; 
	
	// \partial f_{11} / \partial y taken at x = 7, y = 13
	std::cout<< der_funcs::partial_derivative(*f11,x,0.1,2,2) <<"\n";

	// \partial f_{12} / \partial x taken at x = 7, y = 13
	std::cout<< der_funcs::partial_derivative(*f12,x,0.1,2,1) <<"\n";

	// \partial f_{12} / \partial y taken at x = 7, y = 13
	std::cout<< der_funcs::partial_derivative(*f12,x,0.1,2,2) <<"\n";

	x = new double[3+1]; 
	x[1] = -7.0; x[2] = 3.0; x[3] = -5.0; 

	// \partial g_{31} / \partial x taken at x = -7, y = 3, z = -5
	std::cout<< der_funcs::partial_derivative(*g31,x,0.1,3,1) <<"\n";

	// \partial g_{32} / \partial y taken at x = -7, y = 3, z = -5
	std::cout<< der_funcs::partial_derivative(*g32,x,0.1,3,2) <<"\n";

	// \partial g_{33} / \partial z taken at x = -7, y = 3, z = -5
	std::cout<< der_funcs::partial_derivative(*g33,x,0.1,3,3) <<"\n"<<"\n";
}

void examples::test_coord_func_output()
{
	/* Test the outputs of the coordinate functions */

	find_root_multi_D test1(f11,f12); 

	double *x = new(double [2+1]); 

	x[1] = 7.0; x[2] = 3.0; 

	test1.fill_F(x, true);

	find_root_multi_D test2(g11, g12, g13); 

	x = new(double [3+1]); 

	x[1] = -7.0; x[2] = 3.0; x[3] = -5.0; 

	test2.fill_F(x, true);

	delete[] x; 
}

void examples::test_jacobi_output()
{
	// run a test to ensure that the Jacobi matrix is being computed correctly

	find_root_multi_D test2(g21, g22, g23); 

	double *x = new(double [3+1]); 

	x[1] = 0.0; x[2] = 3.0; x[3] = 5.0; 

	test2.fill_F(x, true);

	test2.jacobi_mat(x, true);
}

void examples::test_2D1()
{
	// Find the roots of the function F(X) = 0, with F = { f_{11}(x, y), f_{12}(x, y) }^{T}, X = {x, y}

	std::cout<<"/*****************************************************/\n";
	std::cout<<"2D Function 1\n"; 

	find_root_multi_D func2D1(f11,f12); 

	double *x = new(double [2 + 1]); 

	x[1] = 1.0; x[2] = 1.0; 

	func2D1.newton_raphson_search(x,1.0e-6);
	func2D1.broyden_search(x,1.0e-6);
	func2D1.gradient_search(x,1.0e-6);
	func2D1.gradient_broyden_search(x,1.0e-6);

	delete[] x; 
}

void examples::test_2D2()
{
	// Find the roots of the function F(X) = 0, with F = { f_{21}(x, y), f_{22}(x, y) }^{T}, X = {x, y}

	std::cout<<"\n/*****************************************************/\n";
	std::cout<<"2D Function 2\n";

	find_root_multi_D func2D2(f21,f22); 

	double *x = new(double [2 + 1]); 
	
	x[1] = 2.0; x[2] = 2.0; 

	func2D2.newton_raphson_search(x,1.0e-6);
	func2D2.broyden_search(x,1.0e-6);
	func2D2.gradient_search(x,1.0e-6);
	func2D2.gradient_broyden_search(x,1.0e-6);

	delete[] x; 
}

void examples::test_2D3()
{
	// Find the roots of the function F(X) = 0, with F = { f_{31}(x, y), f_{32}(x, y) }^{T}, X = {x, y}

	std::cout<<"\n/*****************************************************/\n";
	std::cout<<"2D Function 3\n";

	double *x = new(double [2 + 1]); 

	find_root_multi_D func2D3(f31,f32); 

	x[1] = 2.0; x[2] = 2.0; 

	func2D3.newton_raphson_search(x,1.0e-6);
	func2D3.broyden_search(x,1.0e-6);
	func2D3.gradient_search(x,1.0e-6);
	func2D3.gradient_broyden_search(x,1.0e-6);

	delete[] x; 
}

void examples::test_2D4()
{
	// Find the roots of the function F(X) = 0, with F = { f_{41}(x, y), f_{42}(x, y) }^{T}, X = {x, y}

	std::cout<<"\n/*****************************************************/\n";
	std::cout<<"2D Function 4\n";

	double *x = new(double [2 + 1]); 

	find_root_multi_D func2D4(f41,f42); 

	x[1] = 1.0; x[2] = 2.0; 

	func2D4.newton_raphson_search(x,1.0e-6);
	func2D4.broyden_search(x,1.0e-6);
	func2D4.gradient_search(x,1.0e-6);
	func2D4.gradient_broyden_search(x,1.0e-6);

	delete[] x; 
}

void examples::test_3D1()
{
	// Find the roots of the function F(X) = 0, with F = { g_{31}(x, y, z), g_{32}(x, y, z), g_{33}(x, y, z) }^{T}, X = {x, y, z}

	std::cout<<"\n/*****************************************************/\n";
	std::cout<<"3D Function 1\n";

	find_root_multi_D func3D1(g11,g12,g13); 

	double *X = new(double [3 + 1]); 

	X[1] = 0.0; X[2] = 0.0; X[3] = 0.0; 

	func3D1.newton_raphson_search(X,1.0e-6);
	func3D1.broyden_search(X,1.0e-6);
	func3D1.gradient_search(X,1.0e-6);
	func3D1.gradient_broyden_search(X,1.0e-6);

	delete[] X; 
}

void examples::test_3D2()
{
	// Find the roots of the function F(X) = 0, with F = { g_{21}(x, y, z), g_{22}(x, y, z), g_{23}(x, y, z) }^{T}, X = {x, y, z}

	std::cout<<"\n/*****************************************************/\n";
	std::cout<<"3D Function 2\n";

	find_root_multi_D func3D1(g21,g22,g23); 

	double *X = new(double [3 + 1]); 

	X[1] = 0.0; X[2] = 0.0; X[3] = 0.0; 

	func3D1.newton_raphson_search(X,1.0e-6);
	func3D1.broyden_search(X,1.0e-6);
	func3D1.gradient_search(X,1.0e-6);
	func3D1.gradient_broyden_search(X,1.0e-6);

	delete[] X; 
}

void examples::test_3D3()
{
	// Find the roots of the function F(X) = 0, with F = { g_{31}(x, y, z), g_{32}(x, y, z), g_{33}(x, y, z) }^{T}, X = {x, y, z}

	std::cout<<"\n/*****************************************************/\n";
	std::cout<<"3D Function 3\n";

	find_root_multi_D func3D1(g31,g32,g33); 

	double *X = new(double [3 + 1]); 

	X[1] = 0.0; X[2] = 0.0; X[3] = 0.0; 

	func3D1.newton_raphson_search(X,1.0e-6);
	func3D1.broyden_search(X,1.0e-6);
	func3D1.gradient_search(X,1.0e-6);
	func3D1.gradient_broyden_search(X,1.0e-6);

	delete[] X; 
}

void examples::test_3D4()
{
	// Find the roots of the function F(X) = 0, with F = { g_{41}(x, y, z), g_{42}(x, y, z), g_{43}(x, y, z) }^{T}, X = {x, y, z}

	std::cout<<"\n/*****************************************************/\n";
	std::cout<<"3D Function 4\n";

	find_root_multi_D func3D1(g41,g42,g43); 

	double *X = new(double [3 + 1]); 

	X[1] = 1.0; X[2] = 1.0; X[3] = 1.0; 

	func3D1.newton_raphson_search(X,1.0e-6);
	func3D1.broyden_search(X,1.0e-6);
	func3D1.gradient_search(X,1.0e-6);
	func3D1.gradient_broyden_search(X,1.0e-6);

	delete[] X; 
}