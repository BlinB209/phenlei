#pragma once
#ifdef __cplusplus
extern "C"
{
#endif

//! In mgrid, the variables types are defined as:
//!   typedef int idxtype;
//!   typedef short idxtype;
typedef int idxtype;
typedef double realtype;


void MGridGen(int, idxtype *, realtype *, realtype *, idxtype *, realtype *,
              int, int, int *, int *, int *, idxtype *);

#ifdef __cplusplus
}
#endif