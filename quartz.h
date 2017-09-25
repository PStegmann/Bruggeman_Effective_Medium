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

#pragma once

#include <complex>
#include <vector>
#include <fstream>		 /* istringstram */
#include <iostream>
#include <sstream>		 /* stringstream */
#include <string>		   /* string line */
#include <stdlib.h> 	 /* atof */

namespace minerals
{
  class quartz  // CLASS DECLARATION
  {
  private:
    bool debug;
    std::vector<double> database_l;
    std::vector<double> database_n;
    std::vector<double> database_lk;
    std::vector<double> database_k;
  public:

    quartz();

    ~quartz();

    std::complex<double> operator()(double);
  };
}
