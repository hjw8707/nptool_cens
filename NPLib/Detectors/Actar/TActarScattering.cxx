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
 *  This class hold Actar Scattring events properties                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TActarScattering.h"

ClassImp(TActarScattering)

TActarScattering::TActarScattering():
Theta(0),
Phi(0),
Range(0),
Charge(0),
PhiShort(0),
ThetaShort(0),
Voxels(0),
VoxelsInTrigger(0),
Energy(0)
{}