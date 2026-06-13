/* Copyright (C)
 * 2019 - Hu Ren, rh890127a@163.com
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */
/**
 * @file utilityUsingCpp.hpp
 * @brief Declare some c++ standard feathers
 * @author Hu Ren, rh890127a@163.com
 * @version v0.1
 * @date 2019-08-16
 */

#ifndef UTILITY_USINGCPP_HPP
#define UTILITY_USINGCPP_HPP

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#if __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1700)
#include <unordered_map>
#include <unordered_set>
#else
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#endif

#include "stdlib.h"

namespace UTILITY {
/**
 * IO stream
 */
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;
using std::ostream;
using std::streambuf;
using std::stringstream;

/**
 * std container
 */
using std::map;
using std::multimap;
using std::multiset;
using std::pair;
using std::set;
using std::unordered_map;
using std::unordered_set;
using std::vector;

/**
 * std string
 */
using std::string;
/**
* @brief using this to replace std::to_string that add extra ZEROs
         in float number conversion ( noticed in icc 2018 );
         ZEROs before effective figure will be neglected
*/
template <class U>
inline string to_string(const U data) {
  stringstream sstrm;
  // preserve enouph precision for float value
  // commonly 16 for double
  sstrm.precision(16);
  sstrm << data;
  return sstrm.str();
}

/**
 * @brief specilization for float
 */
template <>
inline string to_string<float>(const float data) {
  stringstream sstrm;
  // preserve enouph precision for float value
  // commonly 7 for float
  sstrm.precision(7);
  sstrm << data;
  return sstrm.str();
}

/**
 * transform string to other variable types
 */
#if __cplusplus < 201103L

// To preserve presicion, atof is actually returning double
//#define stof(str) strtof((str).c_str(), NULL)
#define stof(str) atof((str).c_str())
#define stod(str) strtod((str).c_str(), NULL)
#define stold(str) strtold((str).c_str(), NULL)
#define stoi(str) atoi((str).c_str())
#define stol(str) atol((str).c_str())
#define stoll(str) atoll((str).c_str())

#else

// To preserve presicion, stof is actually returning double
using std::stod;
using std::stof;
using std::stoi;
using std::stol;
using std::stold;
using std::stoll;

#endif

}  // namespace UTILITY

#endif  // HSF_USINGCPP_HPP
