#ifndef SWMACRO_H
#define SWMACRO_H

#include "utilities.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BLOCKNUM64K 64
#define EPS 1e-6

#ifdef DEBUG
#define LOG(format, ...) \
  printf("File: " __FILE__ ",Line: %05d: " format "\n", __LINE__, ##__VA_ARGS__)
#else
#define LOG(format, ...)
#endif

#define SLAVE_FUNC(funcname) slave_##funcname
#define ALIGNED(addr) (((((unsigned long)(addr)-1) >> 5) + 1) << 5)
#define ArraySize 57344

// typedef int label;
// typedef int label32;
// typedef long label64;

// typedef double scalar;
// typedef float scalar32;
// typedef double scalar64;

#ifdef __cplusplus
}
#endif

#endif
