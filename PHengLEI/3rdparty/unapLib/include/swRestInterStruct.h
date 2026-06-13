#ifndef SWRESTINTERSTRUCT_H
#define SWRESTINTERSTRUCT_H

#include "swMacro.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  const label *mapPtr;    // restriction map
  const scalar *fPtr;     // fine data
  scalar *cPtr;           // coarse data
  label **localStartEnd;  // local range of fine and coarse data in each slave
                          // core
  label slaveCycles;  // slave cycles (how many cycles of slave cores to store
                      // the entire array)
} restStruct;

typedef struct {
  label *mapPtr;          // interpolation map
  label *offsetMapPtr;    // offset map
  scalar *fPtr;           // fine data
  const scalar *cPtr;     // coarse data
  label **localStartEnd;  // local range of fine and coarse data in each slave
                          // core
  label slaveCycles;  // slave cycles (how many cycles of slave cores to store
                      // the entire array)
} interStruct;

typedef struct {
  const label *mapPtr;    // restriction map
  const scalar *fPtr;     // fine upper data
  scalar *cUPtr;          // coarse upper data
  scalar *cDPtr;          // coarse diagonal data
  label **localStartEnd;  // local range of fine and coarse data in each slave
                          // core
  label slaveCycles;  // slave cycles (how many cycles of slave cores to store
                      // the entire array)
} aggMatrixUpperStruct;

void restrictData_host(restStruct *);

void SLAVE_FUNC(restrictData_slave)(restStruct *);

void interpolateData_host(interStruct *);

void SLAVE_FUNC(interpolateData_slave)(interStruct *);

#ifdef __cplusplus
}
#endif

#endif
