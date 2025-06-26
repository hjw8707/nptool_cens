#ifndef NPTIMESTAMP_h
#define NPTIMESTAMP_h
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author : Hugo Jacob hjacob@ijclab.in2p3.fr                       *
 *                                                                           *
 * Creation Date   : May 2024                                                *
 * Last update     : May 2024                                                *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Resources useful for TS treatment                                        *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>

// NPL
#include "NPInputParser.h"
using namespace NPL;

// ROOT header
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

namespace NPL {

  class TimeStamp {

  public:
    TimeStamp();
    ~TimeStamp(){cout << "test dest TS" << endl;};

   public: // Various Method
    void ReadConfigurationFile();
    // First argument is TS name, second one is TS value
    void AddTimeStamp(std::string TSname, unsigned long long TS){TSMap[TSname] = TS;};
    // Clears the map of TS
    void ClearTimeStamps(){ this->TSMap.clear();};
    //  
    void FindPrompts();
    //
    bool MatchTS(map<std::string,unsigned long long>::iterator TS);
    bool MatchTS(std::string TS);

    std::map<std::string, unsigned long long> GetGoodTS(){return GoodTS;};
   private:
    // Parameters to initialize by reading the config
    // First Parameter is mean of the prompt, second argument is width accepted for the prompt of the TS diff
    std::map<std::string, std::pair<unsigned long long, unsigned int>> TSConditions;
    
    // Map of all timestamps added
    std::map<std::string, unsigned long long> TSMap;
    // Initialization of 1 TS that will be the reference for all other TS
    std::string TSReference;
    // Vector of all the good TS after FindPrompts()
    std::map<std::string, unsigned long long> GoodTS;
    
    int fVerboseLevel;

    ClassDef(TimeStamp, 1)
  };
} // namespace NPL
#endif
