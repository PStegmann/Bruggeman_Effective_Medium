/*
*
*	KK_analysis.cpp
*
* Description:
* ============
*
*	2016-09-07:
*	Test code for the QUAGS algorithm in the GSL
*	library used to solve the integral
*	\int_0^1 x^{-1/2} log(x) dx = -4
*	by numerical quadrature, which has an algebraic
*	logarithmic singularity at the origin.
*
*	2016-09-08:
*	The modified version computes the improper integral
*	\int_{\inf}^{\inf} 1 / ( 1 + x^2 ) dx = pi
*	using the qagi subroutine.
*   
*   2017-03-31:
*   Code for computing the real part of a refractive
*   index spectrum from the imaginary part using
*   trapezoidal quadrature.
*
*	Author:				    Date:		    Revisions:
*	=======				    =====		    ==========
*	Patrick Stegmann		2016-09-07	Original Code
*	Patrick Stegmann		2016-09-08	Makefile repaired
*									   	qags rewritten to qagi
* 	Patrick Stegmann  		2016-10-13  rewriting to qawc for KK analysis.
* 	Patrick Stegmann  		2016-10-19  Applying trapezoidal rule for integration
*                               		instead (i.e. admitting defeat).
*   Patrick Stegmann        2017-03-31  Reads the imaginary part of the
*                                       refractive index from file "refrindex_imag.txt"
*                                       and computes the Kramers-Kronig transform
*                                       using the trapezoidal rule.
*
*/

/*
 * Copyright © 2016, 2017 Patrick Stegmann
 *
 * This file is part of Bruggeman_Effective_Medium.
 *
 * Bruggeman_Effective_Medium is free software:
 * you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*** MACROS ***/
#define LUDOLF 3.14159265359
#define LIGHTSPEED 299792458
#define FACTOR 1.

/*** HEADER ***/
#include <stdio.h>
#include <vector>
#include <math.h>
#include <fstream>		 			/* istringstram */
#include <iostream>
#include <sstream>		 			/* stringstream */
#include <string>		   			/* string line 	*/
//#include <gsl/gsl_integration.h>	/* quadrature 	*/
/*** Boost headers ***/
#include <boost/filesystem.hpp> // filesystem access
#include <boost/assign/std/vector.hpp>

/* linear interpolation */
template<typename T1, typename T2>
inline T1 lin_interp(T1 Q11, T1 Q22, T2 x, T2 x1, T2 x2)
{
  return (Q11 +
		 (x-x1)
		*(Q22-Q11)
		/(x2-x1));
} /* lin_interp end */

/*** INTEGRAL KERNEL  ***/

inline long double f ( double x, void *params)
{
	/*
	*
	* input: - wavelength
	*		 - pointer to wavelengths and corresponding
	*          Imaginary part of refractive index
	*
	* output: Kramers-Kronig integral kernel at a specific
	*         wavelength for quadrature
	*
	*/

	// cast *void pointer to *double
	double *intern = (double *)params;
	double c = *(intern+1500); 	// pass singularity location as last value
								// of parameter pointer
	double t = x;
	t = 1./t;
	c = 1./c;
	//std::cout << c << std::endl;

	// Test output of pointer
	/*
	for (size_t ii = 0; ii < 1001; ii++)
	{
		//beta[ii] = *intern;
		std::cout << *(double *)intern << std::endl;
		intern++;
	}
	*/
	
		// use wavelength to interpolate corresponding Imaginary refractive index.
	double wavelength[1500];
	for (size_t ii = 0; ii < 1500; ii++) 
	{
		wavelength[ii] = ii*0.1;
		//wavelength[ii] = 1./(1+wavelength[ii]);
	}

	double k=0.;
	double kc=0;

	if (t<wavelength[0]) 
	{
		k = lin_interp(*intern,0.,t,wavelength[0],0.);
		//k = *intern;
	} 
	else 
	{
		if (t>=wavelength[1499]) 
		{
			/*k = lin_interp(*(intern+998),
										*(intern+999),
										t,
										wavelength[998],
										wavelength[999]);*/
			//k = *(intern+999);
      		k = *(intern+1499)*exp(-(t-wavelength[1499])/10000);
		} 
		else 
		{
			for (size_t jj = 0; jj < 1500; jj++) 
			{
				if (t >= wavelength[jj] && t<wavelength[jj+1]) 
				{
					k = lin_interp(*(intern+jj),
												*(intern+jj+1),
												t,
												wavelength[jj],
												wavelength[jj+1]);
					break;
				}
			}
		}
	}

	if (c<wavelength[0])
	{
		kc = lin_interp(*intern,0.,c,wavelength[0],0.);
		//kc = *intern;
	}
	else
	{
		if (c>=wavelength[1499]) 
		{
				/*kc = lin_interp(*(intern+998),
										*(intern+999),
										c,
										wavelength[998],
										wavelength[999]);*/
			//kc = *(intern+999);
      		kc = *(intern+1499)*exp(-(c-wavelength[1499])/10000);
      	}
	 	else
	  	{
			for (size_t jj = 0; jj < 1500; jj++)
			{
				if (c >= wavelength[jj] && c<wavelength[jj+1])
				{
					kc = lin_interp(*(intern+jj),
									*(intern+jj+1),
									c,
									wavelength[jj],
									wavelength[jj+1]);
					break;
				}
			}
		}
	}

	t = FACTOR*t;
	c = FACTOR*c;
	t = 1./t;
	c = 1./c;
	// actual function evaluation
	//double f = (t*t*k-c*t*kc)/(t*t-c*c);
	long double f = (t*t*k-c*t*kc)/(t*t-c*c);
	//long double f = k/(t*(c*c-t*t));
	return f;
} /* end of integral kernel f */


/*******************|
*					|
*	MAIN FUNCTION 	|
*					|
********************/

int main (void)
{
	bool debug = 0;

	std::ifstream size_source;
    std::vector<double> size_vector;
    size_source.open("./size_bin.dat", std::ios_base::in);  // open data
    if(!size_source)
    {
      std::cerr << "Size input file not found!\n";
    }
    else
    {
      size_t counter = 1;
      for (std::string line; getline(size_source,line); ) // read stream line by line
      {
        if(counter > 1) // stream also reads bullshit lines at the beginning and end.
        {
          std::istringstream in(line);    // make a stream for the line itself

          double x; // input buffer for database
          in >> x; // actual reading from file
          if ( debug ) std::cout << x << " " << x << " "; // output for debbuging
          size_vector.push_back(x);
        }
        counter += 1;
      }
      if ( debug ) std::cout << "number of lines of the aerosol size table: " << counter << std::endl;
    }
    size_source.close();


    for (std::vector<double>::iterator it = size_vector.begin() ; it != size_vector.end(); ++it)
    {
    	std::string tempFolder("./Index_");
	    std::ostringstream ss;
	    ss << *it;
	    tempFolder += ss.str();
	    boost::filesystem::path outputPath ( tempFolder );
	    boost::system::error_code returnedError;
	    boost::filesystem::path workingDirectory ( boost::filesystem::current_path() );
	    boost::filesystem::create_directories( outputPath);
	    boost::filesystem::current_path(outputPath);
		/****	READ Imaginary refractive index k database ****/
		std::vector<double> database_k, database_r;
		std::ifstream source_k; 									// build a read-Stream
		source_k.open("refrindex_imag.txt", std::ios_base::in);	// open data
		std::ifstream source_r; 									// build a read-Stream
		source_r.open("refrindex_real.txt", std::ios_base::in);	// open data
		std::ofstream output_file;
		output_file.open("refrindex_real_kk.txt");
		size_t counter = 0;
		if(!source_k)
		{
			std::cerr << "Refractive index (imaginary part) data file not found!\n";
		}
		else
		{
			for (std::string line; getline(source_k,line); ) 	// read stream line by line
			{
				if(counter >= 0) 								// stream also reads bullshit lines at the beginning and end.
				{
					std::istringstream in(line);				// make a stream for the line itself

					double x; 		// input buffer for database (size 3 for 3 columns to read)
					in >> x;		// actual reading from file, column 1
					if ( debug ) std::cout << x << std::endl;	// output for debbuging
					database_k.push_back(x); // odd index -> imaginary part of refractive index
				}
				counter += 1;
			}
			if ( 1 ) std::cout << "number of lines of the imaginary database: " << counter << std::endl;
		}
		if(!source_r)
		{
			std::cerr << "Refractive index (real part) data file not found!\n";
		}
		else
		{
			for (std::string line; getline(source_r,line); ) 	// read stream line by line
			{
				if(counter >= 0) 							// stream also reads bullshit lines at the beginning and end.
				{
					std::istringstream in(line);			// make a stream for the line itself

					double x; 		// input buffer for database (size 3 for 3 columns to read)
					in >> x;		// actual reading from file, column 1
					if ( debug ) std::cout << x << std::endl;	// output for debbuging
					database_r.push_back(x); // odd index -> imaginary part of refractive index
				}
			}
		}
		source_k.close();
		source_r.close();
	
		std::cout << boost::filesystem::current_path() << std::endl;
	  	double *database_array = &database_k[0]; // double pointer is linked to vector.
		long double result=0., error;
		
		/*** UV correction ***/
		double nu = 1.;
		database_k.push_back(nu); // pushback wavelength to parameter array
		long double realpart=0., singularity = 0.;
		size_t intervals = 100000;
		double h=(log(6000)-log(0.0000000000000000000001))/intervals;
		
		/*** TRAPEZOIDAL QUADRATURE ***/
		for (size_t ii = 0; ii < intervals; ii++) {
			double lower = exp(log(0.0000000000000000000001)+ii*h);
			double upper = exp(log(0.0000000000000000000001)+(ii+1)*h);
			// SINGULARITY INTERVAL CHECK
			if (nu < upper && nu >= lower)
			{
				double a,b,t,s,k;
				s = nu - lower;
				t = upper - nu;
				k = database_k[10];
				// exact solution at the singularity
				singularity = b*(t+s)+0.5*(a+k-b*nu)*log((2*nu+t)/(2*nu-s));
				//std::cout << "singularity " << singularity << std::endl;
				continue;	
			}
			else
			{
				// trapezoidal rule
				realpart += 0.5*(log(upper)-log(lower))*(f(lower,&database_array[0])+f(upper,&database_array[0]));
			}
		}
		result = 1. + 2.f/LUDOLF*realpart + singularity;
		double UV_correction = result-database_r[10];
		database_k.pop_back();
		realpart = 0.;
		singularity = 0.;

	  	/*** SPECTRAL SAMPLING LOOP ***/
	  	for (size_t kk = 1; kk < counter; kk++)
	  	{
	  		if ( debug ) std::cout << counter << " counter " << std::endl;
			double nu = 1./(kk*0.1);
			if (1) std::cout << kk << " Singularity: " << nu << std::endl;
			database_k.push_back(nu); // pushback wavelength to parameter array
			long double realpart=0., singularity = 0.;
			size_t intervals = 100000;
			double h=(log(6000)-log(0.0000000000000000000001))/intervals;
			 
			/*** TRAPEZOIDAL QUADRATURE ***/
			for (size_t ii = 0; ii < intervals; ii++) {
				double lower = exp(log(0.0000000000000000000001)+ii*h);
				double upper = exp(log(0.0000000000000000000001)+(ii+1)*h);
				if ( debug ) std::cout << intervals << lower << upper << std::endl;
				// SINGULARITY INTERVAL CHECK
				if (nu < upper && nu >= lower)
				{
					double a,b,t,s,k;
					s = nu - lower;
					t = upper - nu;
					k = database_k[kk];
					// exact solution at the singularity
					singularity = b*(t+s)+0.5*(a+k-b*nu)*log((2*nu+t)/(2*nu-s));
					if ( debug ) std::cout << "singularity " << singularity << std::endl;
					continue;	
				}
				else
				{
					// trapezoidal rule
					realpart += 0.5*(log(upper)-log(lower))*(f(lower,&database_array[0])+f(upper,&database_array[0]));
					if ( debug ) std::cout << realpart << " " << kk << " " << ii << std::endl;
				}
			}
			result = 1. + 2.f/LUDOLF*realpart + singularity - UV_correction;
			output_file << result << std::endl;
			database_k.pop_back();
		}

		output_file.close();
		boost::filesystem::current_path(workingDirectory);
  	}
  	return EXIT_SUCCESS;
}
