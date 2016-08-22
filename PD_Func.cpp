#ifndef ATTACH_H
#include "Attach.h"
#endif

double der_funcs::partial_derivative(double (*the_func)(double *x, int n), double *x, double dx, int n, int dirn)
{
	// Use Richardson extrapolation to estimate the partial derivative of a
	// function at the point x
	// this implementation is based on the dfridr algorithm given in
	// NRinC by Press et al. 
	// this implementation uses less function evaluations than the 
	// one provided in the lectures 

	// the_func is a vector valued real function evaluated at position vector x = [x_{1}, x_{2}, ...., x_{n}]

	try{
	
		// argument dirn tells the function the axis along which the derivative is to be computed
		// i.e. dirn = 1 => compute \partial f / \partial x_{1}, dirn = 2 => compute \partial f / \partial x_{2}, etc...
		// so dirn must by necessity be in the range 1 <= dirn <= n

		bool c1 = n > 0 ? true : false; 
		bool c2 = dirn >=0 && dirn <= n ? true : false; 
		bool c3 = dx > 1.0e-6 ? true : false; // ensure that the initial dx is positive to avoid division by zero

		if(c1 && c2 && c3){

			const int ntab = 10; 

			//http://www.cplusplus.com/reference/new/nothrow/
			//http://www.cplusplus.com/reference/new/bad_alloc/			 

			double *xphh = new(std::nothrow)(double [n+1]); 
			double *xmhh = new(std::nothrow)(double [n+1]); 

			double a[ntab][ntab]; // keep the function values in a table

			bool c4 = xphh != nullptr ? true : false; 
			bool c5 = xmhh != nullptr ? true : false; 
			bool c6 = a != nullptr ? true : false; 

			if(c4 && c5 && c6){

				const double con = 1.4, con2 = template_funcs::DSQR(con); 
				const double big = std::numeric_limits<double>::max(); 
				const double safe = 2.0; 

				int i, j, k; 
				double err, errt, fac, hh, ans;

				hh = dx;

				// Need a vector version of x +/- hh
				for(k=1; k<=n; k++){
					xphh[k] = x[k]; xmhh[k] = x[k]; 
				}
				xphh[dirn] += hh; xmhh[dirn] -= hh;

				a[0][0] = (the_func(xphh,n)-the_func(xmhh,n)) / (2.0*hh); // first approximation to f'(x)

				err = big; 

				for(i=1; i<ntab; i++){

					hh /= con; 

					// Need a vector version of x +/- hh
					for(k=1; k<=n; k++){
						xphh[k] = x[k]; xmhh[k] = x[k]; 
					}
					xphh[dirn] += hh; xmhh[dirn] -= hh;
			
					a[0][i] = (the_func(xphh, n)-the_func(xmhh, n)) / (2.0*hh); // approximation to f'(x) with smaller step size
			
					fac = con2; 
			
					// extrapolate the derivative to higher orders without extra function evaluations
					for(j=1; j<=i; j++){

						a[j][i] = (a[j-1][i]*fac-a[j-1][i-1]) / (fac-1.0); 
				
						fac = con2*fac; 
				
						errt = std::max(fabs(a[j][i]-a[j-1][i]), fabs(a[j][i]-a[j-1][i-1])); 
				
						// compute the new error with the error from the previous step
						if(errt <= err){
							err = errt; 
							ans = a[j][i]; // update the derivative estimate
						}
					}

					// if error has increased significantly stop
					if(fabs(a[i][i]-a[i-1][i-1]) >= safe*err){
						break; 
					}

				}

				//std::cout<<"The value of the derivative at x = "<<x<<" is "<<ans<<"\n"; 

				delete[] xphh;
				delete[] xmhh; 

				return ans; 				
			}
			else{
				std::string reason = "Error: double der_funcs::partial_derivative(double (*the_func)(double *x, int n), double *x, double dx, int n, int dirn)\n"; 
				if(!c4) reason += "No memory available for xphh\n"; 
				if(!c5) reason += "No memory available for xmhh\n"; 
				if(!c6) reason += "No memory available for a\n";

				throw std::runtime_error(reason);
			}	
		}
		else{
			std::string reason = "Error: double der_funcs::partial_derivative(double (*the_func)(double *x, int n), double *x, double dx, int n, int dirn)\n"; 
			if(!c1) reason += "Size of position array is input as negative n = " + template_funcs::toString(n) + "\n"; 
			if(!c2) reason += "Direction along which derivative should be computed is not clear dirn = " + template_funcs::toString(dirn) + "\n";
			if(!c3) reason += "Initial spacing value is negative, or too small, dx = " + template_funcs::toString(dx, 9) + "\n";

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