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
//! @file      Geo_OversetGridTopo_Struct.h
//! @brief     It defines the topology information of the overlap grid.
//! @author    Zhang Jian, He Xin.

#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{
//! @brief Geo_OversetGridTopo_Struct class
class Geo_OversetGridTopo_Struct
{
private:
    int numberOfCores;
    int zoneStartPointLabel;
    int *iLinkPointLabel, *jLinkPointLabel, *kLinkPointLabel;
    int *iBlank;
    int zoneStartCenterLabel;
    int *hingedPointContainer;
    RDouble *xCoreContainer, *yCoreContainer, *zCoreContainer;

public:
    Geo_OversetGridTopo_Struct();
    ~Geo_OversetGridTopo_Struct();
public:

    int GetNumberOfCores() const;

    int GetZoneStartPointLabel() const;

    int * GetILinkPointLabel() const;

    int * GetJLinkPointLabel() const;

    int * GetKLinkPointLabel() const;

    int * GetCellTypeContainer() const;

    int * GetHingedPointContainer() const;

    RDouble * GetXCoreContainer() const;

    RDouble * GetYCoreContainer() const;

    RDouble * GetZCoreContainer() const;

    void SetNumberOfCores(int numberOfCores);

    void SetCellContainer(int *iBlank) { this -> iBlank = iBlank; }

    void SetZoneStartPointLabel(int zoneStartPointLabel);

    void SetZoneStartCenterLabel(int zoneStartCenterLabel);

    void CreateLinkPointStructure(int nTotalNode);

    void CreateCellContainer(int nTotalCell);

    void CreateCoreParameterContainer(int numberOfCores);
};

#include "Geo_OversetGridTopo_Struct.hxx"
}