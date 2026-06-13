/**
 * @file: utilities.h
 * @author: Liu Hongbin
 * @brief: some useful utilities
 * @date:   2019-10-08 15:12:44
 * @last Modified by:   lenovo
 * @last Modified time: 2019-12-03 16:10:48
 */
#ifndef UTILITY_UTILITIES_H
#define UTILITY_UTILITIES_H

#include "mpi.h"
#include "rcmf.h"
#include "utilityBasicFunction.h"
#include "utilityExceptions.h"
#include "utilityInterfaces.h"
#include "utilityTimers_c.h"
#include "utilityType.h"

#ifdef __cplusplus

#include <algorithm>
#include <iostream>
#include <vector>

#include "utilityCommunicationManager.hpp"
#include "utilityCommunicator.hpp"
#include "utilityContainer.hpp"
#include "utilityDummyOStream.hpp"
#include "utilityMpiWrapper.hpp"
#include "utilityMultiOStream.hpp"
#include "utilityOStream.hpp"
#include "utilityTimers.hpp"
#include "utilityUsingCpp.hpp"
#include "utilityWriteToFile.hpp"

using namespace UTILITY;
#endif

// #include "utilInterfaces.h"
// #include "exceptions.h"

// } // end namespace HSF

#endif
