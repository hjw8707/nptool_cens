#ifndef TMUST2PHYSICSREADER_H
#define TMUST2PHYSICSREADER_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr      *
 *                                                                           *
 * Creation Date  : 2023                                                     *
 * Last update    : 2023
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  TTreader class for Must2 Physics                                         *
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
#include "TMust2Data.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TTree.h"
#include "NPVTreeReader.h"


using namespace std;

class TMust2PhysicsReader : public NPL::VTreeReader {
public:
  TMust2PhysicsReader();
  ~TMust2PhysicsReader(){};
  
public:
  void r_SetTreeReader(TTreeReader* TreeReader);
private:
  TTreeReader *dummy = new TTreeReader();
public:
  TTreeReaderValue<TMust2Data>* r_ReaderEventData = new TTreeReaderValue<TMust2Data>(*dummy,"");

public:
  ClassDef(TMust2PhysicsReader,0);

};


#endif
