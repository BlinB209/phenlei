#ifndef DRAND48_H
#define DRAND48_H
#if defined(_WIN32) || defined(WIN32)

#include <stdlib.h>
#include <time.h>
#define mmm 0x100000000LL
#define ccc 0xB16
#define aaa 0x5DEECE66DLL

double drand48();
void srand48(unsigned int i);

#endif
#endif