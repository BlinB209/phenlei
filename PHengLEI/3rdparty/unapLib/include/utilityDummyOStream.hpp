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
 * @file DummyOStream.hpp
 * @brief an output stream that will act as "do nothing"
 * @author Hu Ren, rh890127a@163.com
 * @version v0.1
 * @date 2019-08-18
 */

#ifndef UTILITY_DUMMYOSTREAM_HPP
#define UTILITY_DUMMYOSTREAM_HPP

#include "utilityOStream.hpp"

// using namespace UTILITY;

namespace UTILITY {
class DummyOStream : public OStream {
 public:
  //--------------------------------------------------------------
  // construct & deconstruct
  //--------------------------------------------------------------

  // construct empty
  DummyOStream() : OStream() {}

  // deconstruct
  virtual ~DummyOStream() {}

  //--------------------------------------------------------------
  // streaming operator
  //--------------------------------------------------------------
  virtual DummyOStream &operator<<(char chrt) {
    // do nothing
    return *this;
  }

  virtual DummyOStream &operator<<(string str) {
    // do nothing
    return *this;
  }

  virtual DummyOStream &operator<<(int64_t val) {
    // do nothing
    return *this;
  }

  virtual DummyOStream &operator<<(int32_t val) {
    // do nothing
    return *this;
  }

  virtual DummyOStream &operator<<(unsigned long val) {
    // do nothing
    return *this;
  }

  virtual DummyOStream &operator<<(unsigned int val) {
    // do nothing
    return *this;
  }

  virtual DummyOStream &operator<<(double val) {
    // do nothing
    return *this;
  }

  virtual DummyOStream &operator<<(float val) {
    // do nothing
    return *this;
  }

  /**
   * @brief operator<<, interface accept OsOp type parameters
   * @param[in] opt, represent parameter like "ENDL" and "FLUSH".
   * @return
   */
  virtual DummyOStream &operator<<(OsOp opt) {
    // do nothing
    return *this;
  }
};

}  // namespace UTILITY

#endif  // HSD_DUMMYOSTREAM_HPP
