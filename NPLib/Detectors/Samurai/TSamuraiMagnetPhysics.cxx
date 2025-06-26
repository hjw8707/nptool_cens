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
 *  This class hold SamuraiMagnet treated data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TSamuraiMagnetPhysics.h"

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

ClassImp(TSamuraiMagnetPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TSamuraiMagnetPhysics::TSamuraiMagnetPhysics(){
  }

///////////////////////////////////////////////////////////////////////////
void TSamuraiMagnetPhysics::Clear(){
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TSamuraiMagnetPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SAMURAIMagnet");
  // nothing to do
}

///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TSamuraiMagnetPhysics::GetSpectra() {
  map< string , TH1*> empty;
  return empty;
} 

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSamuraiMagnetPhysics::Construct(){
  return (NPL::VDetector*) new TSamuraiMagnetPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_samuraiMagnet{
    public:
      proxy_samuraiMagnet(){
        NPL::DetectorFactory::getInstance()->AddToken("Samurai","Samurai");
        NPL::DetectorFactory::getInstance()->AddDetector("Samurai",TSamuraiMagnetPhysics::Construct);
      }
  };

  proxy_samuraiMagnet p_samuraiMagnet;
}

