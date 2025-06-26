/*****************************************************************************
 * Copyright (C) 2009-2017   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Damien Thisse  contact address: damien.thisse@cea.fr     *
 *                                                                           *
 * Creation Date  : November 2024                                            *
 * Last update    : November 2024                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Actar Beam properties                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TActarBeam.h"

ClassImp(TActarBeam)

TActarBeam::TActarBeam():
    f_EntryPoint(0,0,0),
    f_ExitPoint(0,0,0),
    f_Voxels(0)
{}

void TActarBeam::SetPoints(TVector3 & EntryPoint, TVector3 & ExitPoint){
    f_EntryPoint = EntryPoint;
    f_ExitPoint = ExitPoint;
    f_Theta = (f_ExitPoint - f_EntryPoint).Angle(TVector3(1,0,0));
}

void TActarBeam::SetVoxels(int Voxels){
    f_Voxels = Voxels;
}

void TActarBeam::SetIsScatteredBeam(bool IsScatteredBeam){
    f_IsScatteredBeam = IsScatteredBeam;
}