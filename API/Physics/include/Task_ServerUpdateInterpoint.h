//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Task_ServerUpdateInterpoint.h
//! @brief     Communication pattern of update interpoint data.
//!            It is not recommended for the learners in the first stage.
//! @author    Wan Yunbo.
//! @time      2019.12.9.

#pragma once
#include "Task.h"
namespace PHSPACE
{
//! @brief TaskServerUpdateInterpoint define the task of updating data at interpoint
class LIB_EXPORT TaskServerUpdateInterpoint : public Task
{
public:
    //! Construct function
    TaskServerUpdateInterpoint();
    //! Destructor, free memory
    ~TaskServerUpdateInterpoint();

    //! run the non blocking task for updating the interpoint.
    //! @param[in] actkey         the actkey strore the information of the interpoint.
    void MainTaskNonBlocking(ActionKey *actkey);

    //! run the blocking task for updating the interpoint.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void MainTaskBlocking(ActionKey *actkey);

private:
    //! prepare for the task for update the interpoint.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void PreTask(ActionKey *actkey);

    //! post solve for the task for update the interpoint.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void PostTask(ActionKey *actkey);

    //! main solve for the task for update the interpoint.
    //! @param[in] actkey       the actkey store the information of the interpoint.
    void MainTask(ActionKey *actkey);

    //! communicate the interpoint using MPI.
    //! @param[in] actkey       the actkey store the information of the interpoint.
    //! @param[in] sendBuffer       a data container containing the data for sending.
    //! @param[in] receivedBuffer       a data container containing the data for receiving.
    //! @param[in] iZone       the number of sending zone.
    //! @param[in] jZone       the number of receiving zone.
    //! @param[in] requestContainer       a vector of MPI_Request.
    void InterpointNonBlocking(ActionKey *actkey, DataContainer *sendBuffer, DataContainer *receivedBuffer, int iZone, int jZone, vector<PH_Request> &requestContainer);

    //! get the interpoint using MPI.
    //! @param[in] actkey       the actkey store the information of the interpoint.
    //! @param[in] iZone       the number of sending zone.
    //! @param[in] jZone       the number of receiving zone.
    //! @param[out] nlen       the size of the sending or receiving data.
    CharVecSizeType InterpointGetLength(ActionKey *actkey, int iZone, int jZone);
};

}