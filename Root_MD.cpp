#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the class used to find the roots of multidimensional functions
// R. Sheehan 25 - 3 - 2013

find_root_multi_D::find_root_multi_D()
{
	// Default constructor

	F = nullptr; 
	dx = nullptr; 
	xold = nullptr; 
	xnew = nullptr; 
	J = nullptr; 
	func_vec = nullptr; 

	maxit = ndim = 0; 
}

find_root_multi_D::find_root_multi_D(double (*func1)(double *, int ), double (*func2)(double *, int ))
{
	// Constructor
	try{

		ndim = 2; 

		maxit = 100; 

		// Declare the necessary vectors and matrices
		F = new (double [ndim+1]); 
		dx = new (double [ndim+1]); 
		xold = new (double [ndim+1]);  
		xnew = new (double [ndim+1]); 

		J = new (double *[ndim+1]); 
		for(int i=1; i<=ndim; i++){
			J[i] = new (double [ndim+1]); 
		}

		func_vec = new (coordfunc [ndim+1]); 
		func_vec[1] = func1; 
		func_vec[2] = func2;

	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: find_root_multi_D::find_root_multi_D(double (*func1)(double *, int ), double (*func2)(double *, int ))\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

find_root_multi_D::find_root_multi_D(double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int ))
{
	// Constructor
	try{
		ndim = 3; 

		maxit = 100; 

		// Declare the necessary vectors and matrices
		F = new (double [ndim +1]); 
		dx = new (double [ndim +1]); 
		xold = new (double [ndim+1]);  
		xnew = new (double [ndim+1]);  

		J = new (double *[ndim+1]); 
		for(int i=1; i<=ndim; i++){
			J[i] = new (double [ndim+1]); 
		}

		func_vec = new (coordfunc [ndim+1]); 
		func_vec[1] = func1; 
		func_vec[2] = func2; 
		func_vec[3] = func3; 
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: find_root_multi_D::find_root_multi_D(double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int ))\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

// Methods
double find_root_multi_D::sum_sqr(double *x)
{
	// Sum the squares of the coordinate functions at the point x
	// This is used in the gradient search algorithm

	try{

		bool c1 = x != nullptr ? true : false; 
		bool c2 = func_vec != nullptr ? true : false; 
		bool c3 = ndim > 1 ? true : false; 

		if(c1 && c2 && c3){

			double t,s; 

			t = s = 0.0; 

			for(int i=1; i<=ndim; i++){
				t = template_funcs::DSQR( func_vec[i](x, ndim) ); 
				s += t; 
			}

			return s; 
		
		}
		else{
			std::string reason = "Error: double find_root_multi_D::sum_sqr(double *x)\n";
			if(!c1) reason += "vector x has not been allocated correctly\n"; 
			if(!c2) reason += "vector func_vec has not been allocated correctly\n";
			if(!c3) reason += "num. dimensions less than 2: ndim = " + template_funcs::toString(ndim) + "\n"; 
			
			throw std::invalid_argument(reason); 
		}
		
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void find_root_multi_D::fill_F(double *x, bool loud)
{
	// Assign the value of the vector function F at the point x

	try{

		bool c1 = x != nullptr ? true : false; 
		bool c2 = func_vec != nullptr ? true : false; 
		bool c3 = ndim > 1 ? true : false; 

		if(c1 && c2 && c3){

			for(int i=1; i<=ndim; i++){
				F[i] = func_vec[i](x,ndim); 
			}

			if(loud){
				std::cout<<"\n";
				for(int i=1; i<=ndim; i++){
					std::cout<<F[i]<<"\n";
				}
				std::cout<<"\n";
			}
		}
		else{
			std::string reason = "Error: void find_root_multi_D::fill_F(double *x)\n";
			if(!c1) reason += "vector x has not been allocated correctly\n"; 
			if(!c2) reason += "vector func_vec has not been allocated correctly\n"; 
			if(!c3) reason += "num. dimensions less than 2: ndim = " + template_funcs::toString(ndim) + "\n"; 
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void find_root_multi_D::jacobi_mat(double *pos, bool loud)
{
	// Return the jacobi matrix for the mapping F at the point x

	try{

		bool c1 = pos != nullptr ? true : false; 
		bool c2 = func_vec != nullptr ? true : false; 
		bool c3 = ndim > 1 ? true : false; 
		bool c4 = J != nullptr ? true : false; 

		if(c1 && c2 && c3 && c4){

			for(int i=1; i<=ndim; i++){
				for(int j=1; j<=ndim; j++){
					J[i][j] = der_funcs::partial_derivative( func_vec[i], pos, 0.1, ndim, j ); 
				}
			}

			if(loud){
				for(int i=1; i<=ndim; i++){
					for(int j=1; j<=ndim; j++)
						std::cout<<J[i][j]<<" "; 
					std::cout<<"\n";
				}
			}
		
		}
		else{
			std::string reason = "Error: void find_root_multi_D::jacobi_mat(double *pos, bool loud)\n";
			if(!c1) reason += "vector x has not been allocated correctly\n"; 
			if(!c2) reason += "vector func_vec has not been allocated correctly\n";
			if(!c4) reason += "array J has not been allocated correctly\n";
			if(!c3) reason += "num. dimensions less than 2: ndim = " + template_funcs::toString(ndim) + "\n"; 
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void find_root_multi_D::sherman_morrison_inverse(double **Aold, double **Anew, double *dxold, double *dF, bool loud)
{
	// Compute an approximate inverse of a matrix by the Sherman-Morrison formula

	try{

		bool c1 = Aold != nullptr ? true : false; 
		bool c2 = Anew != nullptr ? true : false; 
		bool c3 = ndim > 1 ? true : false; 
		bool c4 = dxold != nullptr ? true : false;
		bool c5 = dF != nullptr ? true : false;

		if(c1 && c2 && c3 && c4 && c5){

			double *u = new(std::nothrow)(double [ndim+1]); 
			double *v = new(std::nothrow)(double [ndim+1]);

			bool c6 = u != nullptr ? true : false;
			bool c7 = v != nullptr ? true : false;

			if(c6 && c7){

				double lambda, t1, t2; 

				// Define the vectors for the perturbation
				lambda = 0.0; 

				for(int i=1; i<=ndim; i++){
	
					t1=0.0; t2=0.0; 
		
					for(int j=1; j<=ndim; j++){

						t1 += Aold[i][j]*dF[j]; // Aold * dF
			
						t2 += dx[j]*Aold[j][i]; // dx^{T}*Aold
					}
		
					u[i] = dx[i] - t1; // u = dx - Aold*dF
		
					v[i] = t2; // v = dx^{T}*Aold
		
					lambda += dx[i]*t1; // dx * Aold * dF
				}

				if(loud){
					std::cout<<"\nlambda = "<<lambda<<"\n";

					std::cout<<"\nu = "; 
					for(int j=1; j<=ndim; j++){
						std::cout<<u[j]<<"  ";
					}

					std::cout<<"\nv = "; 
					for(int j=1; j<=ndim; j++){
						std::cout<<v[j]<<"  ";
					}
				}

				// Define the inverse approximation from the outer product
				for(int i=1; i<=ndim; i++){
					for(int j=1; j<=ndim; j++){
						Anew[i][j] = Aold[i][j] + ( ( u[i]*v[j] ) / lambda ); 
					}
				}

				delete[] u;
				delete[] v;		
			
			}
			else{
				std::string reason = "Error: void find_root_multi_D::sherman_morrison_inverse(double **Aold, double **Anew, double *dxold, double *dF, bool loud)\n"; 
				
				if(!c6) reason += "No memory available for u\n";
				if(!c7) reason += "No memory available for v\n";

				throw std::runtime_error(reason);
			}
		}
		else{
			std::string reason = "Error: void find_root_multi_D::fill_F(double *x)\n";
			if(!c1) reason += "array Aold has not been allocated correctly\n"; 
			if(!c2) reason += "array Anew has not been allocated correctly\n";
			if(!c4) reason += "vector dxold has not been allocated correctly\n";
			if(!c4) reason += "vector df has not been allocated correctly\n";
			if(!c3) reason += "num. dimensions less than 2: ndim = " + template_funcs::toString(ndim) + "\n"; 
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
	catch(std::runtime_error &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void find_root_multi_D::newton_raphson_search(double *xinit, double toler, bool loud)
{
	// Find the roots of F(x) = 0 using the Newton-Raphson search in multi-dim

	try{

		bool c1 = xinit != nullptr ? true : false; 
		bool c2 = toler > 1.0e-16 ? true :false; 

		if(c1 && c2){

			int iter = 1;
			bool cgt = false; 
			double t1,t2; 

			xold = xinit; 

			while(iter < maxit){
		
				// Fill the Jacobi matrix
				jacobi_mat(xold);
		
				// Evaluate the function at the current position
				fill_F(xold);

				for(int i=1; i<=ndim; i++){
					F[i]*=-1.0; 
				}

				// Find the update vector
				dx = lin_slv::solve_system(J,F,ndim);

				// Compute the error
				t1=0.0; 
				for(int i=1; i<=ndim; i++){
					t2=dx[i]; 
					if(fabs(t2)>t1){
						t1 = t2; 
					}
				}

				// Update the position of the root
				for(int i=1; i<=ndim; i++){
					xnew[i] = xold[i] + dx[i]; 
				}

				if(loud){
					std::cout<<"Iteration "<<iter<<"\n";
					for(int i=1; i<=ndim; i++){
						std::cout<<xnew[i]<<"\n";
					}
					std::cout<<"\n";
				}

				// Check for convergence
				if(fabs(t1) < toler){
					std::cout<<"Newton-Raphson Search has converged to a root in "<<iter<<" iterations\n"; 
					for(int i=1; i<=ndim; i++){
						std::cout<<xnew[i]<<" , "<<func_vec[i](xnew,ndim)<<"\n";
					}
					std::cout<<"\n";
					cgt = true; 
					break; 
				}

				xold = xnew; 

				iter++; 

			}

			if(cgt == false){
				std::cout<<"Newton-Raphson Search failed to converge to a root in "<<maxit<<" iterations\n"; 
				for(int i=1; i<=ndim; i++){
					std::cout<<xnew[i]<<" , "<<func_vec[i](xnew,ndim)<<"\n";
				}
				std::cout<<"\n";
			}
		
		}
		else{
			std::string reason = "Error: void find_root_multi_D::newton_raphson_search(double *xinit, double toler)\n";
			if(!c1) reason += "vector xinit has not been allocated correctly\n"; 
			if(!c2) reason += "tolerance = " + template_funcs::toString(toler, 9) + " is too small\n";
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void find_root_multi_D::broyden_search(double *xinit, double toler, bool loud)
{
	// Find the roots of F(x) = 0 using Broyden's search in multi-dim

	try{

		bool c1 = xinit != nullptr ? true : false; 
		bool c2 = toler > 1.0e-16 ? true :false; 

		if(c1 && c2){		 

			double *Fold = new(std::nothrow)(double [ndim+1]); 
			double *dF = new(std::nothrow)(double [ndim+1]);  
			double **Ainv; 
			double **Ainvold = new(std::nothrow)(double *[ndim+1]); 

			for(int j=1; j<=ndim; j++){
				Ainvold[j] = new(std::nothrow)(double [ndim+1]);  
			}

			bool c3 = Fold != nullptr ? true: false; 
			bool c4 = dF != nullptr ? true: false; 
			bool c5 = Ainvold != nullptr ? true: false; 

			if(c3 && c4 && c5){

				int iter = 1;
				bool cgt = false; 
				double t1,t2;

				xold = xinit; 

				while(iter < maxit){
		
					if(iter == 1){
						// Fill the Jacobi matrix and compute its inverse
						jacobi_mat(xold);

						Ainv = lin_slv::find_inverse(J,ndim); 

						if(loud){
							std::cout<<"\n";
							for(int i=1; i<=ndim; i++){
								for(int j=1; j<=ndim; j++)
									std::cout<<Ainv[i][j]<<" ";
								std::cout<<"\n";
							}
						}

						// Evaluate the function at the current position
						fill_F(xold);
					}
					else{
						// Apply the Sherman-Morrison formula for computing the inverse
			
						// Store F_{old} and A^{-1}_{old}
						for(int i=1; i<=ndim; i++){
							Fold[i] = F[i]; 
							for(int j=1; j<=ndim; j++){
								Ainvold[i][j] = Ainv[i][j]; 
							}
						}
			
						fill_F(xold);
			
						for(int i=1; i<=ndim; i++){
							dF[i] = F[i]-Fold[i]; 
						}

						//Ainv = inverse by sherman morrison
						sherman_morrison_inverse(Ainvold,Ainv,dx,dF);

						if(loud){
							std::cout<<"\n";
							for(int i=1; i<=ndim; i++){
								for(int j=1; j<=ndim; j++)
									std::cout<<Ainv[i][j]<<" ";
								std::cout<<"\n";
							}
						}

					}

					// Find the update vector
					// dx = Ainv*F
					for(int i=1; i<=ndim; i++){
						t1=0.0; 
						for(int j=1;j<=ndim; j++){
							t1 += Ainv[i][j]*F[j]; 
						}
						dx[i] = -1.0*t1; 
					}

					// Update the position of the root
					for(int i=1; i<=ndim; i++){
						xnew[i] = xold[i] + dx[i]; 
					}

					if(loud){
						std::cout<<"\nIteration "<<iter<<"\n";
						for(int i=1; i<=ndim; i++){
							std::cout<<xnew[i]<<"\n";
						}
						std::cout<<"\n";
					}

					// Compute the error
					t1=0.0; 
					for(int i=1; i<=ndim; i++){
						t2=dx[i]; 
						if(fabs(t2)>t1){
							t1 = t2; 
						}
					}

					// Check for convergence
					if(fabs(t1) < toler){
						std::cout<<"Broyden Search has converged to a root in "<<iter<<" iterations\n"; 
						for(int i=1; i<=ndim; i++){
							std::cout<<xnew[i]<<" , "<<func_vec[i](xnew,ndim)<<"\n";
						}
						std::cout<<"\n";
						cgt = true; 
						break; 
					}

					xold = xnew; 

					iter++; 

				}

				delete[] Fold;
				delete[] dF;
				delete[] Ainv; 
				delete[] Ainvold; 

				if(cgt == false){
					std::cout<<"Broyden Search failed to converge to a root in "<<maxit<<" iterations\n"; 
					for(int i=1; i<=ndim; i++){
						std::cout<<xnew[i]<<" , "<<func_vec[i](xnew,ndim)<<"\n";
					}
					std::cout<<"\n";
				}				
			}
			else{
				std::string reason = "Error: void find_root_multi_D::broyden_search(double *xinit, double toler, bool loud)\n"; 
				if(!c3) reason += "No memory available for Fold\n"; 
				if(!c4) reason += "No memory available for dF\n"; 
				if(!c5) reason += "No memory available for Ainvold\n";

				throw std::runtime_error(reason);
			}
		}
		else{
			std::string reason = "Error: void find_root_multi_D::broyden_search(double *xinit, double toler)\n";
			if(!c1) reason += "vector xinit has not been allocated correctly\n"; 
			if(!c2) reason += "tolerance = " + template_funcs::toString(toler, 9) + " is too small\n";
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
	catch(std::runtime_error &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void find_root_multi_D::gradient_search(double *xinit, double toler, bool loud)
{
	// Find the roots of F(x) = 0 using Broyden's search in multi-dim

	try{

		bool c1 = xinit != nullptr ? true : false; 
		bool c2 = toler > 1.0e-16 ? true :false; 

		if(c1 && c2){

			double *z = new(std::nothrow)(double [ndim+1]);
			double *z2 = new(std::nothrow)(double [ndim+1]);

			bool c3 = z != nullptr ? true : false; 
			bool c4 = z2 != nullptr ? true : false; 

			if(c3 && c4){
				int iter = 1;
				bool cgt = false; 
				double t1,t2,g1,g2,g3,gc,gnew,a2,a3,ac,h1,h2,h3;

				xold = xinit; 

				while(iter < maxit){

					// Compute g1
					g1 = sum_sqr(xold); 

					// Fill the Jacobi matrix
					jacobi_mat(xold);
		
					// Evaluate the function at the current position
					fill_F(xold);

					// Compute the vector z
					for(int j=1; j<=ndim; j++){
						t1=0.0; 
						for(int i=1; i<=ndim; i++){
							t1 += J[i][j]*F[i]; 
						}
						z[j] = 2.0*t1; // z = 2*J^{T}*F
					}

					// Compute the 2-norm of z
					t2=0.0; 
					for(int i=1; i<=ndim; i++){
						t2 += template_funcs::DSQR(z[i]); 
					}
					t2 = sqrt(t2); 

					if(loud){
						std::cout<<"\nz =  "<<"\n";
						for(int i=1; i<=ndim; i++){
							std::cout<<z[i]<<"\n";
						}
						std::cout<<"\n"<<"t2 = "<<t2<<"\n";
					}

					if(t2 > toler){

						// Scale z by its 2-norm
						for(int i=1; i<=ndim; i++){
							z[i] /= t2; 
						}

						// Define the interval over which alpha is sought
						a3 = 1.0; 
			
						for(int i=1; i<=ndim; i++){
							z2[i] = xold[i]-z[i]; 
						}
			
						g3 = sum_sqr(z2); 

						while(g3-g1 > toler){

							a3 *= 0.5; 
				
							for(int i=1; i<=ndim; i++){
								z2[i] = xold[i]-a3*z[i]; 
							}
				
							g3 = sum_sqr(z2); 
				
							if(fabs(a3) < 0.5*toler){
								// Not likely to get better than this
								break;
							}
						}

						a2 = 0.5*a3; 
			
						for(int i=1; i<=ndim; i++){
							z2[i] = xold[i]-a2*z[i]; 
						}
			
						g2 = sum_sqr(z2);

						// Perform the interpolation
						h1 = (g2-g1)/a2; 
						h2 = (g3-g2)/(a3-a2); 
						h3 = (h2-h1)/a3; 
						ac = 0.5*(a2-(h1/h3)); /* Critical value of alpha */

						for(int i=1; i<=ndim; i++){
							z2[i] = xold[i]-ac*z[i]; 
						} 

						gc = sum_sqr(z2);

						if(gc < g3){
							for(int i=1; i<=ndim; i++){
								xnew[i] = z2[i]; 
							} 
						}
						else{
							for(int i=1; i<=ndim; i++){
								xnew[i] = xold[i]-a3*z[i]; 
							} 
						}

						gnew = sum_sqr(xnew);

						if(loud){
							std::cout<<"\nIteration "<<iter<<"\n";
							for(int i=1; i<=ndim; i++){
								std::cout<<xnew[i]<<"\n";
							}
							std::cout<<"\n";
						}

						// Compute the error
						/*for(int i=1; i<=ndim; i++){
							dx[i] = xnew[i] - xold[i]; 
						}

						t1=0.0; 
						for(int i=1; i<=ndim; i++){
							t2=dx[i]; 
							if(fabs(t2)>t1){
								t1 = t2; 
							}
						}*/

						// Check for convergence
						if(fabs(gnew-g1) < toler){
							std::cout<<"Gradient Search has converged to a root in "<<iter<<" iterations\n"; 
							for(int i=1; i<=ndim; i++){
								std::cout<<xnew[i]<<" , "<<func_vec[i](xnew,ndim)<<"\n";
							}
							std::cout<<"\n";
							cgt = true; 
							break; 
						}

						xold = xnew; 

						iter++;	
					} 
					else{
						std::cout<<"Iterations will not continue\n"; 
						break; 
					}

				}

				delete[] z;
				delete[] z2;

				if(cgt == false){
					std::cout<<"Gradient Search failed to converge to a root in "<<maxit<<" iterations\n"; 
					for(int i=1; i<=ndim; i++){
						std::cout<<xnew[i]<<" , "<<func_vec[i](xnew,ndim)<<"\n";
					}
					std::cout<<"\n";
				}
				
			}
			else{
				std::string reason = "Error: void find_root_multi_D::gradient_search(double *xinit, double toler, bool loud)\n"; 
				if(!c3) reason += "No memory available for z\n"; 
				if(!c4) reason += "No memory available for z2\n"; 

				throw std::runtime_error(reason);
			}
		}
		else{
			std::string reason = "Error: void find_root_multi_D::gradient_search(double *xinit, double toler)\n";
			if(!c1) reason += "vector xinit has not been allocated correctly\n"; 
			if(!c2) reason += "tolerance = " + template_funcs::toString(toler, 9) + " is too small\n";
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
	catch(std::runtime_error &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void find_root_multi_D::gradient_broyden_search(double *xinit, double toler, bool loud)
{
	// Start the search using the gradient method
	// Finish the search using the Broyden method

	try{

		bool c1 = xinit != nullptr ? true : false; 
		bool c2 = toler > 1.0e-16 ? true :false; 

		if(c1 && c2){

			// initialise the search with gradient method
			gradient_search(xinit,1.0e-3, loud); 
			
			// refine the search with broyden method
			broyden_search(xnew,toler, loud);		
		}
		else{
			std::string reason = "Error: void find_root_multi_D::gradient_broyden_search(double *xinit, double toler)\n";
			if(!c1) reason += "vector xinit has not been allocated correctly\n"; 
			if(!c2) reason += "tolerance = " + template_funcs::toString(toler, 9) + " is too small\n";
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}
