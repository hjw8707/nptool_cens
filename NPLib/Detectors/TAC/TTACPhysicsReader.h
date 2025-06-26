#ifndef TTACPHYSICSREADER_H
#define TTACPHYSICSREADER_H
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
 *  TTreader class for TAC Physics                                           *
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
#include "TTACData.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TTree.h"
#include "NPVTreeReader.h"


using namespace std;


class TTACPhysicsReader : public NPL::VTreeReader {
public:
  TTACPhysicsReader();
  ~TTACPhysicsReader(){};
  
public:
  void r_SetTreeReader(TTreeReader* TreeReader);
private:
  TTreeReader *dummy = new TTreeReader();
public:
  TTreeReaderValue<TTACData>* r_ReaderEventData = new TTreeReaderValue<TTACData>(*dummy,"");

public:
  ClassDef(TTACPhysicsReader,0);

};


#endif
