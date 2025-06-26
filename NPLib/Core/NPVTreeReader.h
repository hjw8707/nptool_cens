#ifndef VTreeReader_H
#define VTreeReader_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien MATTA  contact address: matta@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   :   June 2009                                             *
 * Last update     :   Nov 2013                                              *
 *---------------------------------------------------------------------------*
 * Decription:  Class VTreeReader (virtual) for NPAnalysis                     *
 * All detector class used in NPAnalysis should derived from this virtual    *
 *  class. Those VTreeReader daughter will deal with geometry, calibration and *
 *  basic data treatment.                                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *   See MUST2 array for exemple of VTreeReader derived class                  *
 *                                                                           *
 *   Adding Fill Spectra Method to fill control Histogramm                   *
 *****************************************************************************/

// ROOT
#include"TTreeReader.h"
//#include"NPVDetector.h"
#include "Rtypes.h"
// STL
#include <map>
#include <string>
#include <vector>

// NPL
#include "NPCore.h"
#include "NPInputParser.h"

namespace NPL {

  class VTreeReader{
  //class VTreeReader : public VDetector {
   public:
    //  Default Constructor and destructor
    VTreeReader();
    virtual ~VTreeReader();
    
    virtual void r_BuildPhysicalEvent();

    //virtual void BuildSimplePhysicalEvent(){};

    virtual void r_ClearEventPhysics();
    virtual void r_ClearEventData();
    virtual void r_SetTreeReader(TTreeReader*){};
    virtual void r_InitializeRootInputRaw(){};
    virtual void r_InitializeRootInputPhysics(){};  
    virtual void r_InitializeRootOutput(){};  
    private:
      ClassDef(VTreeReader,0);
  };

} // namespace NPL

#endif
