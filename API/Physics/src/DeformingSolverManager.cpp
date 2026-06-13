#include "PHMpi.h"
#include "PHIO.h"
#include "Constants.h"
#include "AleManager.h"
#include "Force.h"
#include "Data_Field.h"
#include "Math_BasisFunction.h"
#include "Zone.h"
#include "SixDofManager.h"
#include "DeformingSolverManager.h"

#include <iomanip>
#include <iostream>

namespace PHSPACE
{

DeformingParameter::DeformingParameter()
{
    numberOfParameters = 0;
    //parameters = 0; 
}

DeformingParameter::~DeformingParameter()
{
    //delete [] parameters;
}

DeformingParameter & DeformingParameter::operator = (const DeformingParameter & rhs)
{
    if (this == & rhs) return * this;

    numberOfParameters = rhs.numberOfParameters;
    parameters = rhs.parameters;

    return * this;
}

bool DeformingParameter::operator == (const DeformingParameter & rhs)
{
    if (this == & rhs) return true;

    if (numberOfParameters != rhs.numberOfParameters)
    {
        return false;
    }

    const RDouble eps = 1.0e-7;
    for (int iParameter = 0; iParameter < numberOfParameters; ++ iParameter)
    {
        if (PHSPACE::ABS(parameters[ iParameter ] - rhs.parameters[ iParameter ]) > eps)
        {
            return false;
        }
    }
    return true;
}

void DeformingParameter::SetParameters(PHVectorRDouble1D& parameterVector)
{
    //delete [] parameters;

    numberOfParameters = static_cast<int>(parameterVector.size());
    parameters.resize(numberOfParameters);

    for (int iParameter = 0; iParameter < numberOfParameters; ++ iParameter)
    {
        parameters[ iParameter ] = parameterVector[ iParameter ];
    }
}

DeformingParameter * GetDeformingParameter()
{
    DeformingParameter * tmp = new DeformingParameter;

    return tmp;
}

DeformingSolverManager::DeformingSolverManager()
{

}

DeformingSolverManager::~DeformingSolverManager()
{

}

void DeformingSolverManager::Initialize()
{

}

void DeformingSolverManager::Restart()
{

}

void DeformingSolverManager::DumpRestart()
{

}

void DeformingSolverManager::Run()
{

}

void DeformingSolverManager::Post()
{

}

void DeformingSolverManager::Dump()
{

}

}