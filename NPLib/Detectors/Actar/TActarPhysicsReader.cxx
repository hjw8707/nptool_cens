/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Damien THISSE  contact address: damien.thisse@cea.fr     *
 *                                                                           *
 * Creation Date  : 5th November 2024                                        *
 * Last update    : 5th November 2024                                        *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "TActarPhysicsReader.h"
#include "NPDetectorFactory.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TActarPhysicsReader)

TActarPhysicsReader::TActarPhysicsReader()
{
};

void TActarPhysicsReader::r_SetTreeReader(TTreeReader* TreeReader){
    r_ReaderEventData = new TTreeReaderValue<MEventReduced>(*TreeReader, "data");    
}; 