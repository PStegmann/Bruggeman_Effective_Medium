/*
* get_refractive_index_ice.cpp
*
* Description:
* ============
* Computes the refractive index of Saharan or Asian Dust.
*
* Record of Revisions:
* ====================
*
* Author:       Date:         Modification:
* =======       =====         =============
* P. Stegmann   2016-06-29    Original sceleton code
*                             & makefile
* P. Stegmann	2016-06-30	  get_index Functor
*                             & index file I/O
*							                & matrix class
* P. Stegmann	2016-07-01	  Added wavelength search,
*							                & temperature search
*							                & bilinear interpolation
*							                & Maxwell-Garret for ice/air mix
* P. Stegmann	2016-07-05	  Cleanup and bilinear template f.
* P. Stegmann   2016-09-28    Re-using the code for the calculation of Sand
*                             instead of ice.
* P. Stegmann   2016-09-30    Creating all classes in the namespace minerals.
* P. Stegmann   2016-10-03    Completed mineral classes with file input and
*                             interpolation of real and imaginary refractive
*                             index.
* P. Stegmann   2016-10-05    Adding ODE solver to compute the Bruggeman compo-
*                             sition. Added Bruggeman_solver method,
*                             streaming_observer method, the Bruggeman_RHS
*                             struct Functor, and linked boost/numeric/odeint.hpp
*                             to this end.
* P. Stegmann   2016-10-11    Added list output to file.
*
*/

/*
 * Copyright © 2016 Patrick Stegmann
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

/*** HEADER ***/
#include <fstream>                  /* istringstram */
#include <iostream>
#include <sstream>                  /* stringstream */
#include <string>                   /* string line */
#include <complex>                  /* complex<T> */
#include <stdlib.h>                 /* atof */
#include <valarray>                 /* class matrix */
#include <numeric>                  /* std::accumulate */
/*** Boost headers ***/
#include <boost/numeric/odeint.hpp> /* ODE integrator for Bruggeman solver */
#include <boost/filesystem.hpp> // filesystem access
#include <boost/assign/std/vector.hpp>

// MINERAL database
#include "WeNeedMoreMinerals.h"

using namespace std;

// Definitions:

bool debug = 0;

/*
*
*		Mathematical Matrix class for data storage based on std::valarray<T>
*
*/
template <class element_type>
class matrix
{

private:

    valarray<element_type> m_storage;
    size_t m_stride;
    size_t m_height;

public:

    matrix(size_t width, size_t height): m_storage(width*height), m_stride(width), m_height(height) {  };

    ~matrix(){};

    element_type &operator()(size_t row, size_t column)
    {
    		// Slice approach doesn't work:
        // column major
        // return m_storage[std::slice(column, m_height, m_stride)][row];

        // row major
    		//valarray<element_type> requested_row = m_storage[std::slice(row, m_stride, m_height)];
        //return requested_row[column];

    		// inelegant replacement approach:
    		return m_storage[row*m_stride + column];
    };

};


/*
*
*		Main class
*
*/
class get_index // CLASS DECLARATION & DEFINITION
{

private:

	//matrix<double> database = matrix<double>(23,2001);	// matrix to hold Dr. Yang's database

	/* bilinear interpolation template member */
	template<typename T1, typename T2>
	inline T1 bilinear_interp(T1 Q11, T1 Q12, T1 Q21, T1 Q22, T2 x1, T2 x2, T2 y1, T2 y2, T2 x, T2 y)
	{
		return (Q11*(x2-x)*(y2-y)\
		 			+ Q12*(x2-x)*(y-y1)\
		 			+ Q21*(x-x1)*(y2-y)\
		 			+ Q22*(x-x1)*(y-y1))\
					/((x2-x1)*(y2-y1));
	}; /* bilinear_interp end */

  /* Right-Hand side Functor for Bruggeman ODE integrator */
  struct Bruggeman_RHS
  {
    const valarray<complex<double> > ind;
    valarray<double> composition;
    const double size;

    /* linear interpolation */
    template<typename T1, typename T2>
    inline T1 lin_interp(T1 Q11, T1 Q22, T2 x, T2 x1, T2 x2)
    {
      return (Q11 +
              (x-x1)
              *(Q22-Q11)
              /(x2-x1));
    } /* interp end */

    Bruggeman_RHS(valarray<complex<double> > indices, double particle_size)
    : ind(indices), size(particle_size), composition(7)
    {
      // Mineral Concentrations for Northern Sahara per diameter:
      double diameter[10] = {0.16, 0.35, 0.71, 1.6, 3.5, 7.1, 16, 35, 71, 158};
      double hematite[10] = {0.006, 0.005, 0.007, 0.008, 0.006, 0.007, 0.007, 0.008, 0.006, 0.007};
      double quartz[10] = {0.029, 0.013, 0.098, 0.109, 0.101, 0.111, 0.147, 0.211, 0.446, 0.566};
      double silicate[10] = {0.391, 0.413, 0.560, 0.594, 0.626, 0.64, 0.636, 0.624, 0.491, 0.393};
      double calcite[10] = {0.018, 0.02, 0.085, 0.096, 0.114, 0.114, 0.087, 0.062, 0.019, 0.003};
      double sulphate[10] = {0.515, 0.505, 0.138, 0.073, 0.05, 0.036, 0.027, 0.004, 0.002, 0.};
      double soot[10] = {0., 0.,0., 0., 0., 0., 0., 0., 0.};
      double water[10] = {0.041, 0.044, 0.112, 0.121, 0.102, 0.093, 0.096, 0.091, 0.036, 0.031};

      if (particle_size <= diameter[0])
      {
        // Extrapolation below available data:
        //cout << "too small" << diameter[0] << particle_size << endl;
        
        composition[0] = lin_interp(hematite[1],hematite[0],particle_size,diameter[1],diameter[0]);
        composition[1] = lin_interp(quartz[1],quartz[0],particle_size,diameter[1],diameter[0]);
        composition[2] = lin_interp(silicate[1],silicate[0],particle_size,diameter[1],diameter[0]);
        composition[3] = lin_interp(calcite[1],calcite[0],particle_size,diameter[1],diameter[0]);
        composition[4] = lin_interp(sulphate[1],sulphate[0],particle_size,diameter[1],diameter[0]);
        composition[5] = lin_interp(soot[1],soot[0],particle_size,diameter[1],diameter[0]);
        composition[6] = lin_interp(water[1],water[0],particle_size,diameter[1],diameter[0]);
        /*
        composition[0] = hematite[0];
        composition[1] = quartz[0];
        composition[2] = silicate[0];
        composition[3] = calcite[0];
        composition[4] = sulphate[0];
        composition[5] = soot[0];
        composition[6] = water[0];
        */
      }
      else
      {
        if (particle_size > diameter[9])
        {
          // Extrapolation beyond available data:
          //cout << "too large" << endl;
          
          composition[0] = lin_interp(hematite[8],hematite[9],particle_size,diameter[8],diameter[9]);
          composition[1] = lin_interp(quartz[8],quartz[9],particle_size,diameter[8],diameter[9]);
          composition[2] = lin_interp(silicate[8],silicate[9],particle_size,diameter[8],diameter[9]);
          composition[3] = lin_interp(calcite[8],calcite[9],particle_size,diameter[8],diameter[9]);
          composition[4] = lin_interp(sulphate[8],sulphate[9],particle_size,diameter[8],diameter[9]);
          composition[5] = lin_interp(soot[8],soot[9],particle_size,diameter[8],diameter[9]);
          composition[6] = lin_interp(water[8],water[9],particle_size,diameter[8],diameter[9]);
          /*
          composition[0] = hematite[9];
          composition[1] = quartz[9];
          composition[2] = silicate[9];
          composition[3] = calcite[9];
          composition[4] = sulphate[9];
          composition[5] = soot[9];
          composition[6] = water[9];
          */
        }
        else
        {
          // Interpolation of data
          //cout << "in the middle." << endl;
          for (size_t jj = 0; jj < (sizeof(diameter)/sizeof(*diameter)); jj++) {
            if (particle_size > diameter[jj] && particle_size <= diameter[jj+1]) {
              //cout << jj << endl;
              composition[0] = lin_interp(hematite[jj],hematite[jj+1],particle_size,diameter[jj],diameter[jj+1]);
              composition[1] = lin_interp(quartz[jj],quartz[jj+1],particle_size,diameter[jj],diameter[jj+1]);
              composition[2] = lin_interp(silicate[jj],silicate[jj+1],particle_size,diameter[jj],diameter[jj+1]);
              composition[3] = lin_interp(calcite[jj],calcite[jj+1],particle_size,diameter[jj],diameter[jj+1]);
              composition[4] = lin_interp(sulphate[jj],sulphate[jj+1],particle_size,diameter[jj],diameter[jj+1]);
              composition[5] = lin_interp(soot[jj],soot[jj+1],particle_size,diameter[jj],diameter[jj+1]);
              composition[6] = lin_interp(water[jj],water[jj+1],particle_size,diameter[jj],diameter[jj+1]);
              /*
              cout << composition[0] << endl;
              cout << composition[1] << endl;
              cout << composition[2] << endl;
              cout << composition[3] << endl;
              cout << composition[4] << endl;
              cout << composition[5] << endl;
              cout << composition[6] << endl;
              */
              break;
            }
          }
        }
      }
    }; /* Bruggeman_RHS constructor end */

    /* Bruggeman_RHS Functor overload */
    void operator()(const complex<double> &x, complex<double> &dxdt, double t) const
    {
      valarray<complex<double> > alpha(ind.size());
      alpha[0] = 1.;
      for (size_t ii = 1; ii < ind.size(); ii++)
          alpha[ii] = t;
      complex<double> numerator(0.,0.);
      complex<double> denominator(0.,0.);
      for (size_t kk = 1; kk < 7; kk++)
      {
        numerator += composition[kk]*(ind[kk]-x)/(ind[kk]+2.*x);
      }
      for (size_t kk = 0; kk < 7; kk++)
      {
        denominator += composition[kk]*alpha[kk]*ind[kk]/((ind[kk]+2.*x)*(ind[kk]+2.*x));
      }
      dxdt = 1./3.*numerator/denominator;
    };

  }; /* Bruggeman_RHS end */


  struct streaming_observer
  {
    std::ostream& m_out;

    streaming_observer( std::ostream &out ) : m_out( out ) { }

    template< class State >
    void operator()( const State &x , double t ) const
    {
        m_out << t;
        m_out << "\t" << x.real() << "\t" << x.imag() ;
        m_out << "\n";
    }
  }; /* streaming_observer end */


  complex<double> Bruggeman_solver (valarray<complex<double> > indices, double particle_size)
  {
    using namespace boost::numeric::odeint;
    complex<double> index_mg = indices[0];
    const double dt = 0.005;
    typedef runge_kutta4<complex<double> > stepper_type;
    integrate_const(stepper_type(),
                    Bruggeman_RHS(indices,particle_size),
                    index_mg,
                    0.0,
                    1.0,
                    dt /* streaming_observer(cout) */);
    return index_mg;
  }; /* Bruggeman_solver end */

  /*
  complex<double> NewRa (complex<double> index_mg, valarray<complex<double> > indices, double particle_size)
  {
    for (size_t ii = 0; ii < 10; ii++) {
      complex<double> numerator(0.,0.);
      complex<double> denominator(0.,0.);
      for (size_t kk = 0; kk < 7; kk++)
      {
        numerator += composition[kk]*(indices[kk]-index_mg)/(indices[kk]+2.*index_mg);
        denominator += composition[kk]*indices[kk]/((indices[kk]+2.*index_mg)*(indices[kk]+2.*index_mg));
      }
      index_mg = index_mg + 1./3.
                *numerator
                /denominator;
    }
    return index_mg;
  }
  */

public:

	get_index(){}	// CLASS CONSTRUCTOR

	~get_index(){};	// CLASS DESTRUCTOR

	/*** OPERATOR OVERLOAD for Functor to return the desired refractive index ***/
	complex<double> operator()(double wavelength, double particle_size)
	{
		complex<double> index_mg(1,1);

    minerals::silicate sili;
    minerals::calcite calci;
    minerals::ferrum fer;
    minerals::quartz qua;
    minerals::sulphates sul;
    minerals::soot soo;
    minerals::water wat;

    valarray<complex<double> > indices(7);
    indices[0] = fer(wavelength)*fer(wavelength);
    indices[1] = qua(wavelength)*qua(wavelength);
    indices[2] = sili(wavelength)*sili(wavelength);
    indices[3] = calci(wavelength)*calci(wavelength);
    indices[4] = sul(wavelength)*sul(wavelength);
    indices[5] = soo(wavelength)*soo(wavelength);
    indices[6] = wat(wavelength)*wat(wavelength);
    for (size_t ii  = 0; ii  < 7; ii++) {
      if (0) cout << indices[ii] << endl;
      // Safety against inf or nan data.
      if (indices[ii] != indices[ii])
      {
        indices[ii].real(0.);
        indices[ii].imag(0.);
      }
      if (indices[ii].real() < 0.) {
        indices[ii].real(abs(indices[ii].real()));
      }
      if (indices[ii].imag() < 0.) {
        indices[ii].imag(abs(indices[ii].imag()));
      }
    }

    // Step 1: ODE solver to 'slowly turn on the inhomogeneity'.
    index_mg = Bruggeman_solver(indices, particle_size);
    if (index_mg.real()< 0.) {
      if (debug) cout<< wavelength << endl;
    }

    // Step 2: Newton-Raphson iteration until convergence.
    //index_mg = NewRa(index_mg,indices,particle_size);

		return sqrt(index_mg);
	}; /* operator () overload end */

}; /* class get_index */

/*
*
*
*
*		MAIN function
*
*
*
*/
int main ( int argc, char *argv[] )
{
	// WELCOME MESSAGE
	cout << "WELCOME." << endl;
	cout << " This program computes the refractive index of Saharan or Asian dust." << endl;
	cout << " It is based on wavelength and size dependent composition." << endl;
	cout << " This code was written by Patrick Stegmann." << endl;

	// CHECK COMMANDLINE ARGUMENT SANITY
	if ( argc != 3 ) // argc should be 3 for correct execution (first argument is app name)
	{
		// We print argv[0] assuming it is the program name
		cout<<"usage: "<< argv[0] <<" wavelength[micro-m] size[micro-m]\n";
	}
	else // correct number of arguments
	{
		// Requesting refractive index:
		get_index get_index_now;
		complex<double> index_holder(1,1);
    
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
      if(debug) cout << "Refractive index has been found!" << endl;
      std::string tempFolder("./Index_");
      std::ostringstream ss;
      ss << *it;
      tempFolder += ss.str();
      std::cout << tempFolder << std::endl;
      boost::filesystem::path outputPath ( tempFolder );
      boost::system::error_code returnedError;
      boost::filesystem::path workingDirectory ( boost::filesystem::current_path() );
      boost::filesystem::create_directories( outputPath, returnedError );
      boost::filesystem::current_path(outputPath);
  		// Refractive index output:
  		if (debug) cout << "Compound refractive index: " << index_holder << endl; // output needs to be different from (1,1).
  		ofstream output_file_real;
      ofstream output_file_imag;
  		output_file_real.open ("refrindex_real.txt");
      output_file_imag.open ("refrindex_imag.txt");
      boost::filesystem::current_path(workingDirectory);
      for (size_t iii = 0; iii < 1500; iii++) 
      {
        //index_holder = get_index_now((0.1+0.1*iii),atof(argv[2]));
        index_holder = get_index_now((0.1+0.1*iii),*it);
        output_file_real << index_holder.real() << endl;
        output_file_imag << index_holder.imag() << endl;
        std::cout << iii << std::endl;
      }
  		output_file_real.close();
      output_file_imag.close();
    }
	}
  return EXIT_SUCCESS;
}
