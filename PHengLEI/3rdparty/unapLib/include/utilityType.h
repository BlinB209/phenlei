/*
The MIT License

Copyright (c) 2019 Hanfeng GU <hanfenggu@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * @file: utilityBase.h
 * @author: Hongbin Liu
 * @email:
 * @Date:   2019-10-06 10:07:59
 * @Last Modified by:   Hongbin Liu
 * @Last Modified time: 2019-11-14 09:50:18
 */

#ifndef UTILITY_UTILITYTYPE_HPP
#define UTILITY_UTILITYTYPE_HPP

#include <stdint.h>

/**
 * define long integer type
 */
typedef int64_t label64;
/**
 * define integer type
 */
typedef int32_t label32;
/**
 * define short integer type
 */
typedef int16_t label16;

/**
 * define short float type
 */
typedef float scalar32;
/**
 * define long float type
 */
typedef double scalar64;

// 默认使用整型
/**
 * define integer type according to compiling macro
 */
#if defined(LABEL_INT32)
typedef label32 label;
#elif defined(LABEL_INT64)
typedef label64 label;
#else
typedef label32 label;
#endif

// 默认使用double精度
#if defined(SCALAR_FLOAT32)
typedef scalar32 scalar;
#elif defined(SCALAR_FLOAT64)
typedef scalar64 scalar;
#else
typedef scalar64 scalar;
#endif

#endif
