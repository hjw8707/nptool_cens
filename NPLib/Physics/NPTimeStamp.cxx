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


#include "NPTimeStamp.h"

#include <fstream>
#include <iostream>


// ROOT
#include "TF1.h"

ClassImp(TimeStamp)

TimeStamp::TimeStamp(){

};


// Checking if the TS diff is in the prompt with the TS reference within a certain width
bool TimeStamp::MatchTS(map<std::string,unsigned long long>::iterator TS){
  // Case TS not recorded -> no TS recorded for this event in this channel
  auto TSref = &TSReference;
  if(std::find_if(TSMap.begin(), TSMap.end(), [TSref](std::pair<std::string, unsigned long long> Pair){return Pair.first == *TSref;}) != TSMap.end() ){
    // Checking that the channel has been declared in the config
    if(std::find_if(TSConditions.begin(), TSConditions.end(), [TS](std::pair<std::string, std::pair<unsigned long long,unsigned int>> Pair){return Pair.first == (*TS).first;}) != TSConditions.end()
      || (*TS).first == *TSref ){
      // Particular case of TS = 0 -> unusual behavior, cout an error message, possible TS problem
      if(TSMap[(*TS).first] > 0){

        // TS reference always correlated with itself
        if((*TS).first == TSReference)
          return true;
        // Checking if the event is in the prompt with TS reference
        else
          return abs((long long)(TSMap[(*TS).first] -TSMap[TSReference] - TSConditions[(*TS).first].first)) < TSConditions[(*TS).first].second;
      }
    
      // Very dubious case
      else{
        std::cout << "WARNING, TS = 0 recorded for channel " << (*TS).first << std::endl;
        return false;
      }
    }
    else{
      std::cout << "Channel " << (*TS).first << " was not initialized in TSconfig" << std::endl;
      return false;
    }
  }
  // Case of TS == 0 -> assume no event recorded for this channel
  else
    return false;
    // return abs(  (long long) ((*TS).second-TSMap[TSReference] - TSConditions[(*TS).first].first)) < TSConditions[(*TS).first].second;
};

// Checking if the TS diff is in the prompt with the TS reference within a certain width
bool TimeStamp::MatchTS(std::string TS){
  // Case TS not recorded -> no TS recorded for this event in this channel
  auto TSref = &TSReference;
  if(std::find_if(TSMap.begin(), TSMap.end(), [TSref](std::pair<std::string, unsigned long long> Pair){return Pair.first == *TSref;}) != TSMap.end() ){
  // if(find (TSMap.begin(), TSMap.end(), TSReference)!= TSMap.end() ){
    // Particular case of TS = 0 -> unusual behavior, cout an error message, possible TS problem
    if(TSMap[TS] > 0){
      
      // TS reference always correlated with itself
      if(TS == TSReference)
        return true;
      // Checking if the event is in the prompt with TS reference
      else
        return abs((long long)(TSMap[TS] -TSMap[TSReference] - TSConditions[TS].first)) < TSConditions[TS].second;
    }
    
    // Very dubious case
    else{
      std::cout << "WARNING, TS = 0 recorded for channel " << TS << std::endl;
      return false;
    }
  }
  // Case of TS == 0 -> assume no event recorded for this channel
  else
    return false;
};

void TimeStamp::FindPrompts(){
  // Checking that the TS reference exists, we assume that at least 1 TS is the ref and is absolutely necessary
  // in the TS treatment
  GoodTS.clear();
  if(MatchTS(TSReference)){
    // CLear the good TS vector and adding the reference TS
    // GoodTS.push_back(TSReference);
    for(auto TS_it = TSMap.begin(); TS_it != TSMap.end(); TS_it++){
      // If prompt with the reference TS, add the TS to the good TS vector
      // else do nothing
      if(MatchTS(TS_it))
        GoodTS[TS_it->first] = TS_it->second;
    }
  }
  else{
    std::cout << "TS matching not done : TS reference not correctly found" << std::endl;
  }
}
void TimeStamp::ReadConfigurationFile() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/TSConfig.dat";

  // open analysis config file
  ifstream TSConfigFile;
  TSConfigFile.open(FileName.c_str());

  if (!TSConfigFile.is_open()) {
    cout << "No TSConfig.dat found: No parameters loaded"  << endl;
    return;
  }
  cout << "Loading user parameter for TS from TSConfig.dat " << endl;

  // read config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!TSConfigFile.eof()) {
    // Pick-up next line
    getline(TSConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigTS";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      TSConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        TSConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="ADD_TS") {
        TSConfigFile >> DataBuffer;
        std::string TSname = DataBuffer;
        TSConfigFile >> DataBuffer;
        long long TSPrompt = atoll(DataBuffer.c_str());
        TSConfigFile >> DataBuffer;
        unsigned int TSWidth = atoi(DataBuffer.c_str());
        TSConditions[TSname] = std::make_pair(TSPrompt,TSWidth);
        cout << whatToDo << " " << TSname << " with prompt value: " << TSPrompt << " and width : " << TSWidth << endl;
      }
      else if (whatToDo=="TS_REF") {
        TSConfigFile >> DataBuffer;
        TSReference = DataBuffer;
        cout << whatToDo << " " <<  DataBuffer << endl;
      }
      else {
        ReadingStatus = false;
      }
    }
  }
}