/******************************************************************************
 * Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
 * HYPRE Project Developers. See the top-level COPYRIGHT file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 ******************************************************************************/

#define HYPRE_RELEASE_NAME    "HYPRE"
#define HYPRE_RELEASE_VERSION "2.24.0"
#define HYPRE_RELEASE_NUMBER  22400
#define HYPRE_RELEASE_DATE    "2022/02/11"
#define HYPRE_RELEASE_TIME    "00:00:00"
#define HYPRE_RELEASE_BUGS    "https://github.com/hypre-space/hypre/issues"

#define HYPRE_DEVELOP_STRING  "v2.24.0-12-gfa43ea82e"
#define HYPRE_DEVELOP_NUMBER   12
#define HYPRE_DEVELOP_BRANCH  "master"

/* Use long long int for HYPRE_BigInt */
/* #undef HYPRE_MIXEDINT */

/* Use long long int for HYPRE_BigInt and HYPRE_Int*/
/* #undef HYPRE_BIGINT */

/* Use single precision values for HYPRE_Real */
/* #undef HYPRE_SINGLE */

/* Use quad precision values for HYPRE_Real */
/* #undef HYPRE_LONG_DOUBLE */

/* Use complex values */
/* #undef HYPRE_COMPLEX */

/* Debug mode */
/* #undef HYPRE_DEBUG */

/* Define to be the max dimension size (must be at least 3) */
#define HYPRE_MAXDIM 3

/* Use persistent communication */
/* #undef HYPRE_USING_PERSISTENT_COMM */

/* Use hopscotch hashing */
/* #undef HYPRE_HOPSCOTCH */

/* Compile without MPI */
/* #undef HYPRE_SEQUENTIAL */

/* Use HYPRE timing routines */
/* #undef HYPRE_TIMING */

/* Use internal BLAS library */
#define HYPRE_USING_HYPRE_BLAS 1

/* Use internal LAPACK library */
#define HYPRE_USING_HYPRE_LAPACK 1

/* Print HYPRE errors */
/* #undef HYPRE_PRINT_ERRORS */

/* Use OpenMP */
#define HYPRE_USING_OPENMP 1

/* Use Caliper instrumentation */
/* #undef HYPRE_USING_CALIPER */

/* Use if executing on device with CUDA */
/* #undef HYPRE_USING_CUDA */

/* Use if executing on device with SYCL */
/* #undef HYPRE_USING_SYCL */

/* Use cuBLAS */
/* #undef HYPRE_USING_CUBLAS */

/* Use CUDA streams */
/* #undef HYPRE_USING_CUDA_STREAMS */

/* Use cuRAND */
/* #undef HYPRE_USING_CURAND */

/* Use cuSPARSE */
/* #undef HYPRE_USING_CUSPARSE */

/* Use device memory pool */
/* #undef HYPRE_USING_DEVICE_POOL */

/* Use unified memory */
/* #undef HYPRE_USING_UNIFIED_MEMORY */

/* Use device memory without UM */
/* #undef HYPRE_USING_DEVICE_MEMORY */

/* Use if executing on device with OpenMP */
/* #undef HYPRE_USING_DEVICE_OPENMP */

/* Use if executing on GPU device */
/* #undef HYPRE_USING_GPU */

/* Use HIP */
/* #undef HYPRE_USING_HIP */

/* Use NVTX */
/* #undef HYPRE_USING_NVTX */

/* Use oneMLK spasre */
/* #undef HYPRE_USING_ONEMKLSPARSE */

/* Use oneMLK blas */
/* #undef HYPRE_USING_ONEMKLBLAS */

/* Use oneMLK rand */
/* #undef HYPRE_USING_ONEMKLRAND */

/* Use SuperLU_Dist */
/* #undef HYPRE_USING_DSUPERLU */

/* Use SuperLU */
/* #undef HAVE_SUPERLU */

/* Use MPI */
#define HYPRE_HAVE_MPI 1

/* #undef HYPRE_HAVE_MPI_COMM_F2C */

/* Define as follows to set the Fortran name mangling scheme:
 * 0 = unspecified
 * 1 = no underscores
 * 2 = one underscore
 * 3 = two underscores
 * 4 = caps, no underscores
 * 5 = one underscore before and after */
#define HYPRE_FMANGLE 0

/* Define as in HYPRE_FMANGLE to set the BLAS name mangling scheme */
#define HYPRE_FMANGLE_BLAS 0

/* Define as in HYPRE_FMANGLE to set the LAPACK name mangling scheme */
#define HYPRE_FMANGLE_LAPACK 0

/* Define to a macro mangling the given C identifier (in lower and upper
 * case), which must not contain underscores, for linking with Fortran. */
#define HYPRE_F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define HYPRE_F77_FUNC_(name,NAME) name ## __

/* Define to 1 if using host memory only */
#define HYPRE_USING_HOST_MEMORY 1
