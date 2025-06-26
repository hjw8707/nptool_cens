/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold GladMagnet treated data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TGladMagnetPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "NPOptionManager.h"
#include "NPDetectorFactory.h"
#include "NPSystemOfUnits.h"
//   ROOT
using namespace NPUNITS;
///////////////////////////////////////////////////////////////////////////

ClassImp(TGladMagnetPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TGladMagnetPhysics::TGladMagnetPhysics(){
  }

///////////////////////////////////////////////////////////////////////////
void TGladMagnetPhysics::Clear(){
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TGladMagnetPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("GLADMagnet");
  // nothing to do
}

///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TGladMagnetPhysics::GetSpectra() {
  map< string , TH1*> empty;
  return empty;
} 

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TGladMagnetPhysics::Construct(){
  return (NPL::VDetector*) new TGladMagnetPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_gladMagnet{
    public:
      proxy_gladMagnet(){
        NPL::DetectorFactory::getInstance()->AddToken("Glad","Glad");
        NPL::DetectorFactory::getInstance()->AddDetector("Glad",TGladMagnetPhysics::Construct);
      }
  };

  proxy_gladMagnet p_gladMagnet;
}

