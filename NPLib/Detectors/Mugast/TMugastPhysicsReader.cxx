/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : febuary 2009                                             *
 * Last update    : July 2021
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Mugast TreeReader                                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TMugastPhysicsReader.h"
#include "NPDetectorFactory.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TMugastPhysicsReader)

TMugastPhysicsReader::TMugastPhysicsReader(){};

void TMugastPhysicsReader::r_SetTreeReader(TTreeReader* TreeReader) {
  r_ReaderEventData = new TTreeReaderValue<TMugastData>(*TreeReader, "Mugast");
};
