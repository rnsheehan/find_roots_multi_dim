#ifndef ATTACH_H
#include "Attach.h"
#endif

void lin_slv::ludcmp(double **a, int *indx, int size, double *d)
{
	// compute the LU decomposition of a matrix a using Crout's algorithm
	// a is written over by this function
	// upon output a contains L and U

	try{

		bool c1 = a != nullptr ? true : false;
		bool c2 = indx != nullptr ? true : false;
		bool c3 = size > 0 ? true : false; 

		if(c1 && c2 && c3){
			
			int i,imax,j,k;
			
			double big,dum,sum,temp;

			static const double TINY=(1.0e-15);
	
			double *vv = new(double [size+1]); // create a vector for storing values for scaling

			*d=1.0; // This keeps track of any row interchanges that are made

			// store scaling values in vv
			for(i=1; i<=size; i++){
				big=0.0;

				for(j=1; j<=size; j++){

					if( (temp=fabs(a[i][j]) ) > big){
						big = temp;
					}

				}
				if(big == 0.0){
					std::cerr<<"Singular matrix in routine LUDCMP\n";
				}

				vv[i] = 1.0/big; 
			}


			for(j=1; j<=size; j++){

				// Start applying Crout's algorithm for computing the LU Decomposition
				for(i=1;i<j;i++){

					sum = a[i][j];
			
					for(k=1; k<i; k++){
						sum -= a[i][k]*a[k][j];
					}
			
					a[i][j] = sum;
				}

				// Find the pivot element
				big=0.0;
				for(i=j; i<=size; i++){

					sum = a[i][j];
			
					for(k=1; k<j; k++){
						sum -= a[i][k]*a[k][j];
					}
			
					a[i][j] = sum;
			
					if((dum = vv[i]*fabs(sum)) >= big){
						big = dum;
						imax = i;
					}
				}

				// Simulate a row swap if required
				if(j!=imax){

					for(k=1; k<=size; k++){
						dum = a[imax][k];
						a[imax][k] = a[j][k];
						a[j][k] = dum;
					}

					*d = -(*d);
					vv[imax] = vv[j];
				}

				indx[j] = imax;
		
				if(a[j][j] == 0.0){
					a[j][j] = TINY;
				}
		
				// Scale a row by the pivot element
				if(j!=size){

					dum = 1.0/(a[j][j]);

					for(i=j+1; i<=size; i++){
						a[i][j] *= dum;
					}

				}
			}

			delete[] vv;	
		}
		else{
			std::string reason = "Error: void lin_slv::ludcmp(double **a, int *indx, int size, double *d)\n"; 
			if(!c1) reason += "Matrix a has not been assigned correctly"; 
			if(!c2) reason += "Array indx has not been assigned correctly";
			if(!c3) reason += "Size of arrays is not positive n = " + template_funcs::toString(size) + "\n";

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	} 
}

void lin_slv::lubksb(double **a, int *indx, double *b, int size)
{
	// solve the system of linear equations a.x = b by LU Decomposition
	// a is assumed to contain its LU decomposition
	// indx contains any row swap information
	// b contains the solution of the system of equations upon output

	try{

		bool c1 = a != nullptr ? true : false;
		bool c2 = indx != nullptr ? true : false;
		bool c3 = size > 0 ? true : false;
		bool c4 = b != nullptr ? true : false; 

		if(c1 && c2 && c3 && c4){

			int i,ii=0,ip,j;
			double sum;

			// Forward substitution to solve L.y = b
			for(i=1; i<=size; i++){

				ip = indx[i];
				sum = b[ip];
				b[ip] = b[i];
		
				if(ii){
					for(j=ii; j<=i-1; j++){
						sum -= a[i][j]*b[j];
					}
				}
				else if(sum){
					ii = i;
				}
				b[i] = sum;
			}

			// Back substitution to solve U.x = y
			for(i=size; i>=1; i--){
				sum = b[i];
				for(j=i+1; j<=size; j++){
					sum -= a[i][j]*b[j];
				}
				b[i] = sum/(a[i][i]); // store the solution in b
			}
		
		}
		else{
			std::string reason = "Error: void lin_slv::lubksb(double **a, int *indx, double *b, int size)\n"; 
			if(!c1) reason += "Matrix a has not been assigned correctly"; 
			if(!c2) reason += "Array indx has not been assigned correctly";
			if(!c3) reason += "Size of arrays is not positive n = " + template_funcs::toString(size) + "\n";

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void lin_slv::ludet(double **a, double *d, int size)
{
	// compute the determinant of a matrix from its LU decomposition

	try{

		bool c1 = a != nullptr ? true : false;
		bool c3 = size > 0 ? true : false;
		bool c4 = d != nullptr ? true : false;

		if(c1 && c3 && c4){

			int j;

			for(j=1; j<=size; j++){
				*d*=a[j][j];
			}

		}
		else{

			std::string reason = "Error: void lin_slv::ludet(double **a, double *d, int size)\n"; 
			if(!c1) reason += "Matrix a has not been assigned correctly"; 
			if(!c4) reason += "Pointer d has not been assigned correctly";
			if(!c3) reason += "Size of arrays is not positive n = " + template_funcs::toString(size) + "\n";

			throw std::invalid_argument(reason); 
			
		}

	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void lin_slv::luinv(double **a, double **y, double *col, int *indx, int size)
{
	// find the inverse of a matrix from its LU Decomposition
	// inverting by solving a.x=col, where col is a unit vector

	try{
	
		bool c1 = a != nullptr ? true : false;
		bool c2 = y != nullptr ? true : false;
		bool c3 = size > 0 ? true : false;
		bool c4 = col != nullptr ? true : false;
		bool c5 = indx != nullptr ? true : false;

		if(c1 && c2 && c3 && c4 && c5){
			
			int i,j;

			for(j=1; j<=size; j++){

				for(i=1; i<=size ;i++){
					col[i]=0.0;
				}

				col[j]=1.0;
		
				lubksb(a,indx,col,size);
		
				for(i=1; i<=size; i++){
					y[i][j] = col[i];
				}
			}

		}
		else{
			std::string reason = "Error: void lin_slv::luinv(double **a, double **y, double *col, int *indx, int size)\n"; 
			if(!c1) reason += "Matrix a has not been assigned correctly"; 
			if(!c2) reason += "Matrix y has not been assigned correctly"; 
			if(!c4) reason += "Array col has not been assigned correctly";
			if(!c5) reason += "Array indx has not been assigned correctly";
			if(!c3) reason += "Size of arrays is not positive n = " + template_funcs::toString(size) + "\n";

			throw std::invalid_argument(reason);
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

double *lin_slv::solve_system(double **a, double *b, int size)
{
	// Sovle the linear system A.x = b using LU Decomposition

	try{

		bool c1 = a != nullptr ? true : false;
		bool c2 = b != nullptr ? true : false;
		bool c3 = size > 0 ? true : false;
	
		if(c1 && c2 && c3){

			// This variable *d is used to keep track of row interchanges inside the lu decomposition function
			double dtmp = 0.0; 
			double *d = &dtmp;  

			//http://www.cplusplus.com/reference/new/nothrow/
			//http://www.cplusplus.com/reference/new/bad_alloc/

			// Create the necessary arrays
			int *indx = new(std::nothrow)(int [size+1]);
			double *sol = new(std::nothrow)(double [size+1]); 
	
			// Declare 2D matrix to hold the data
			double **lua = new(std::nothrow)(double *[size+1]);
	
			for(int i=1; i<=size; i++){
				lua[i] = new(std::nothrow)(double [size+1]); 
			}

			bool c4 = indx != nullptr ?  true : false; 
			bool c5 = sol != nullptr ?  true : false; 
			bool c6 = lua != nullptr ?  true : false; 

			if(c4 && c5 && c6){
				// copy a into lua, then compute its lu decomposition
				// need to keep a copy of a because it is used elsewhere
				for(int i=1; i<=size; i++){
					indx[i] = 0; 
					sol[i] = b[i]; 
					for(int j=1; j<=size; j++){
						lua[i][j] = a[i][j]; 
					}
				}

				// use memcpy as it is more efficient
				// Not sure if this works given how the arrays are declared
				// memcpy requires a contiguous block of memory
				/*memcpy(sol, b, (size+1)*sizeof(double) );
				memcpy(lua, a, (size+1)*(size+1)*sizeof(double) );*/

				// use copy instead of memcpy
				//std::copy(&a[0][0], &a[0][0]+(size+1)*(size+1), &lua[0][0]);

				// compute the lu decomposition of a, this is stored in lua
				ludcmp(lua,indx,size,d);

				// solve the system using the backsubstitution
				lubksb(lua, indx, sol, size); 

				delete[] indx; 
				delete[] lua; 

				return sol;					
			}
			else{
				std::string reason = "Error: double *lin_slv::solve_system(double **a, double *b, int size)\n"; 
				if(!c4) reason += "No memory available for indx\n"; 
				if(!c5) reason += "No memory available for sol\n"; 
				if(!c6) reason += "No memory available for lua\n";

				throw std::runtime_error(reason);
			}
		}
		else{
			std::string reason = "Error: double *lin_slv::solve_system(double **a, double *b, int size)\n"; 
			if(!c1) reason += "Matrix a has not been assigned correctly"; 
			if(!c2) reason += "Array b has not been assigned correctly"; 
			if(!c3) reason += "Size of arrays is not positive n = " + template_funcs::toString(size) + "\n";

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

double **lin_slv::find_inverse(double **a, int size)
{
	// Perform the steps necessary to compute the inverse of the matrix A
	// A copy of a is made
	// The lu decomposition of a is stored in this array
	// the inverse of a is returned by this function

	try{

		bool c1 = a != nullptr ? true : false;
		bool c3 = size > 0 ? true : false;
		
		if(c1 && c3){

			// This variable *d is used to keep track of row interchanges inside the lu decomposition function
			double dtmp = 0.0; 
			double *d = &dtmp;  

			//http://www.cplusplus.com/reference/new/nothrow/
			//http://www.cplusplus.com/reference/new/bad_alloc/

			// Create the necessary arrays
			int *indx = new(std::nothrow)(int [size+1]);
			double *col = new(std::nothrow)(double [size+1]);

			// Declare matrices to hold the data
			double **lua = new(std::nothrow)(double *[size+1]);
			double **ainv = new(std::nothrow)(double *[size+1]);

			for(int i=1; i<=size; i++){
				lua[i] = new(std::nothrow)(double [size+1]); 
				ainv[i] = new(std::nothrow)(double [size+1]); 
			}

			bool c4 = indx != nullptr ?  true : false; 
			bool c5 = col != nullptr ?  true : false; 
			bool c6 = lua != nullptr ?  true : false; 
			bool c7 = ainv != nullptr ?  true : false; 

			if(c4 && c5 && c6 && c7){

				// copy a into lua, then compute its lu decomposition
				for(int i=1; i<=size; i++){
					indx[i] = 0; 
					for(int j=1; j<=size; j++){
						lua[i][j] = a[i][j]; 
					}
				}

				// use copy instead of memcpy
				//std::copy(&a[0][0], &a[0][0]+(size+1)*(size+1), &lua[0][0]);

				// compute the lu decomposition of a, this is stored in lua
				ludcmp(lua,indx,size,d); 

				// compute the inverse of the matrix a using the lu decomposition
				luinv(lua,ainv,col,indx,size);

				delete[] indx;
				delete[] col; 
				delete[] lua;

				return ainv; 		
			
			}
			else{
				std::string reason = "Error: double *lin_slv::solve_system(double **a, double *b, int size)\n"; 
				if(!c4) reason += "No memory available for indx\n"; 
				if(!c5) reason += "No memory available for col\n"; 
				if(!c6) reason += "No memory available for lua\n";
				if(!c7) reason += "No memory available for lua\n";

				throw std::runtime_error(reason);
			}
		}
		else{
			std::string reason = "Error: double **lin_slv::find_inverse(double **a, int size)\n"; 
			if(!c1) reason += "Matrix a has not been assigned correctly"; 
			if(!c3) reason += "Size of arrays is not positive n = " + template_funcs::toString(size) + "\n";

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

//double **lin_slv::sherman_morrison_inverse(double **a, double *u, double *v, int size)
//{
//	// Compute an approximation to the inverse of the matrix (a+u*v), * is an outer product
//	// The inverse is computed via the Sherman-Morrison formula
//
//
//}