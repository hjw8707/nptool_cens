#ifndef TCATSPHYSICSREADER_H
#define TCATSPHYSICSREADER_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
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
 *  This class hold must2 treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <map>
#include <vector>
// NPL
#include "TCATSData.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TTree.h"
#include "NPVTreeReader.h"


using namespace std;

// Forward Declaration

class TCATSPhysicsReader : public NPL::VTreeReader {
public:
  TCATSPhysicsReader();
  ~TCATSPhysicsReader(){};
  
public:
  void r_SetTreeReader(TTreeReader* TreeReader);

private:
  TTreeReader *dummy = new TTreeReader();

public:
  TTreeReaderValue<TCATSData>* r_ReaderEventData = new TTreeReaderValue<TCATSData>(*dummy,"");

public:
  ClassDef(TCATSPhysicsReader,0);

};


#endif
