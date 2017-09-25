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

#include "silicate.h"

minerals::silicate::silicate()
{
  debug = 0;
  /****	READ Silicate refractive index database ****/

  std::ifstream source; 							// build a read-Stream
  source.open("./real/realpart_Silicates.txt", std::ios_base::in);	// open data
  if(!source)
  {
    std::cerr << "Refractive index data file not found!\n";
  }
  else
  {
    size_t counter = 1;
    for (std::string line; getline(source,line); ) // read stream line by line
    {
      if(counter > 1) // stream also reads bullshit lines at the beginning and end.
      {
        std::istringstream in(line);		// make a stream for the line itself

        double x[2]; // input buffer for database
        in >> x[0];	// actual reading from file
        in >> x[1];	// actual reading from file
        if ( debug ) std::cout << x[0] << " " << x[1] << " ";	// output for debbuging
        database_l.push_back(x[0]);
        database_n.push_back(x[1]);
      }
      counter += 1;
    }
    if ( debug ) std::cout << "number of lines of the silicates database: " << counter << std::endl;
  }
  source.close();

  source.open("./imag/imagpart_Silicates.txt", std::ios_base::in);	// open data
  if(!source)
  {
    std::cerr << "Refractive index data file not found!\n";
  }
  else
  {
    size_t counter = 1;
    for (std::string line; getline(source,line); ) // read stream line by line
    {
      if(counter > 1) // stream also reads bullshit lines at the beginning and end.
      {
        std::istringstream in(line);		// make a stream for the line itself

        double x[2]; // input buffer for database
        in >> x[0];	// actual reading from file
        in >> x[1];	// actual reading from file
        if ( debug ) std::cout << x[0] << " " << x[1] << " ";	// output for debbuging
        database_lk.push_back(x[0]);
        database_k.push_back(x[1]);
      }
      counter += 1;
    }
    if ( debug ) std::cout << "number of lines of the silicates database: " << counter << std::endl;
  }
  source.close();

  if (debug) std::cout << "2. Silicates database read from file." << std::endl;
};

minerals::silicate::~silicate(){};

std::complex<double> minerals::silicate::operator()(double wavelength)
{
  // Initialise output container
  std::complex<double> index(1,1);
  if (debug) std::cout << "database size: " << database_l.size() << std::endl;
  // Fetch real part
  if (wavelength < *database_l.begin() ) // Extrapolation below available data
  {
    if (debug) std::cout << "Extrapolation below " << *(database_l.begin()+5) << database_l[*database_l.begin()+5] << std::endl;
    float temp=1.f;
    temp = wavelength*database_n[*database_l.begin()]/database_l[*database_l.begin()];
    /*
    temp = database_n[*database_l.begin()+1] +
          (wavelength-database_l[*database_l.begin()])/(database_l[*database_l.begin()]-database_l[*database_l.begin()+1])
          *(database_n[*database_l.begin()]-database_n[*database_l.begin()+1]);
          */
    index.real(temp);
  }
  else
  {
    if (wavelength >= database_l.back()) // Extrapolation beyond available data
    {
      if (debug) std::cout << "Extrapolation beyond " << database_l[*database_l.end()] << " " << database_l.back() << std::endl;
      float temp=1.f;
      temp = database_n.back()*exp(-(wavelength-database_l.back()/500));
        /*
      temp = database_n.back() +
            (wavelength-database_lk.back())/(database_lk.back()-*(database_lk.end()-2))
            *(database_n.back()-*(database_n.end()-2));
            */
      index.real(temp);
    }
    else // Interpolation in the array of available data.
    {
      if (debug) std::cout << "Interpolation " << std::endl;
      for (ssize_t oo = 0; oo < database_l.size(); oo++) {
        float temp;
        if (wavelength >= database_l[oo] && wavelength < database_l[oo+1]) {
          float temp=1.f;
          temp = database_n[oo] +
                (wavelength-database_l[oo])/(database_l[oo+1]-database_l[oo])
                *(database_n[oo+1]-database_n[oo]);
          index.real(temp);
          break;
        }
      }
    }
  }

  // Fetch imaginary part

  if (wavelength < database_lk[*database_lk.begin()] ) // Extrapolation below available data
  {
    float temp=1.f;
    temp = wavelength*database_k[*database_lk.begin()]/database_lk[*database_lk.begin()];
    /*
    temp = database_k[*database_lk.begin()+1] +
          (wavelength-database_lk[*database_lk.begin()])/(database_lk[*database_lk.begin()]-database_lk[*database_lk.begin()+1])
          *(database_k[*database_lk.begin()]-database_k[*database_lk.begin()+1]);
          */
    index.imag(temp);
  }
  else
  {
    if (wavelength >= database_lk.back()) // Extrapolation beyond available data
    {
      float temp=1.f;
      temp = database_k.back()*exp(-(wavelength-database_lk.back()/500));
      /*
      temp = database_k.back() +
            (wavelength-database_lk.back())/(database_lk.back()-*(database_lk.end()-2))
            *(database_k.back()-*(database_k.end()-2));
            */
      index.imag(temp);
    }
    else // Interpolation in the array of available data.
    {
      if (debug) std::cout << "Interpolation" << std::endl;
      for (ssize_t oo = 0; oo < database_l.size(); oo++) {
        float temp;
        if (wavelength >= database_lk[oo] && wavelength < database_lk[oo+1]) {
          float temp=1.f;
          temp = database_k[oo] +
                (wavelength-database_lk[oo])/(database_lk[oo+1]-database_lk[oo])
                *(database_k[oo+1]-database_k[oo]);
          index.imag(temp);
          break;
        }
      }
    }
  }

  // Output
  return index;
}
