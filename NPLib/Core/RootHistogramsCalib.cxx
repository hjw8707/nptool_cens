/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 21/07/09                                                 *
 * Last update    : 10/05/21                                                 *
 *---------------------------------------------------------------------------*
 * Decription: This class is a singleton class which deals with the ROOT     *
 *             output file and tree both for NPSimulation and NPAnalysis.    *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *   + 03/02/11: Add support for TAsciiFile objects (N. de Sereville)        *
 *   + 10/05/21: Add support for split tree output (A. Matta)                *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "RootHistogramsCalib.h"
#include "NPOptionManager.h"

using namespace std;

RootHistogramsCalib* RootHistogramsCalib::instance = 0;
////////////////////////////////////////////////////////////////////////////////
RootHistogramsCalib* RootHistogramsCalib::getInstance(std::string fileNameBase){
  // A new instance of RootHistogramsCalib is created if it does not exist:
  if (instance == 0) {
    instance = new RootHistogramsCalib(fileNameBase.c_str());
  }

  // The instance of RootHistogramsCalib is returned:
  return instance;
}

////////////////////////////////////////////////////////////////////////////////
void RootHistogramsCalib::Destroy(){ 
  if (instance != 0) {
    delete instance;
    instance = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
RootHistogramsCalib::RootHistogramsCalib(std::string fileNameBase){

  cout << endl << "/////////// ROOT Output files ///////////" << endl;
  cout << "Initializing ouput files " << endl;
  pCurrentDirectory= gDirectory;
  bool calibration=false;

  if(NPOptionManager::getInstance()->IsCalibration()){
    calibration= true;
  }

  // Setup the base name
  if(calibration)
    pBaseName = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  else if(NPOptionManager::getInstance()->IsEfficiency())
    pBaseName = NPOptionManager::getInstance()->GetEfficiencyOutputPath();
  else
    pBaseName="./";
  
  if(fileNameBase.at(0)!='/' && pBaseName.at(pBaseName.length() - 1)!='/')
    pBaseName += "/"+fileNameBase;
  else
    pBaseName+=fileNameBase;


  if (fileNameBase.find("root")==std::string::npos) 
    pBaseName += ".root";

  CreateFile("global");
}

////////////////////////////////////////////////////////////////////////////////
void RootHistogramsCalib::CreateFile(std::string name){
  // Create the tree only if does not exist already
  string file_name=pBaseName;

  cout << " - Creating output file " << file_name.c_str() << endl;
  pRootFile = new TFile(file_name.c_str(), "RECREATE");
  pRootFile->SetCompressionLevel(1);

  gDirectory->cd(pCurrentDirectory->GetPath()); 
  //
  /*
     for(auto it = m_DetectorMap.begin() ; it!=m_DetectorMap.end() ;++it){
     string insertion = "_"+it->first;
     master << filebase << insertion << ".root" << std::endl;
     string filename=path+"/"+filebase+"/"+filebase+insertion+".root";
     auto file = new TFile(filename.c_str(),"RECREATE");
     string treename = "RawTree_"+it->first;
     auto tree = new TTree("RawTree",treename.c_str());
     m_TreeMap[it->first]=tree;
     m_FileMap[it->first]=file;
     tree->SetDirectory(file);
     std::cout << "Splitting tree: " << filename << std::endl;
     it->second->InitBranch(tree);
     }
     master.close();

*/
}
////////////////////////////////////////////////////////////////////////////////
void RootHistogramsCalib::Fill(){
// for(auto it = pRootTrees.begin();it!=pRootTrees.end();it++){
//  it->second->Fill();
// } 
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
RootHistogramsCalib::~RootHistogramsCalib(){ 
  // The data is written to the file and the tree is closed:
  /*if (pRootFiles.size()>0) {
    cout << endl << endl << "Root Output summary" << endl;
    TDirectory* pCurrentDirectory= gDirectory;
    for(auto it = pRootFiles.begin(); it!=pRootFiles.end();it++){
      cout << " - " <<it->first << " tree and file " << endl;
      gDirectory->cd(it->second->GetPath());
      // write TAsciiFile if used
      // EventGenerator
      if (!pEventGenerator->IsEmpty()) pEventGenerator->Write(0,TAsciiFile::kOverwrite);
      // DetectorConfiguration
      if (!pDetectorConfiguration->IsEmpty()) pDetectorConfiguration->Write(0,TAsciiFile::kOverwrite);
      // CalibrationFile
      if (!pCalibrationFile->IsEmpty()) pCalibrationFile->Write(0,TAsciiFile::kOverwrite);
      // RunToTreatFile
      if (!pRunToTreatFile->IsEmpty()) pRunToTreatFile->Write(0,TAsciiFile::kOverwrite);
      // Analysis ConfigFile
      if (!pAnalysisConfigFile->IsEmpty()) pAnalysisConfigFile->Write(0,TAsciiFile::kOverwrite);

      cout << "  -> Number of entries in the " << it->first << " Tree: " << pRootTrees[it->first]->GetEntries() << endl;
      cout << "  -> Number of bites written to file: " << pRootTrees[it->first]->Write(0, TObject::kOverwrite) << endl;
      it->second->Flush();
      it->second->Purge(1);

      gDirectory->cd(pCurrentDirectory->GetPath());
      it->second->Close();
    }
  }

  else {
    cout << "\033[1;31mNo histograms and Tree !\033[0m" << endl;
  }
*/
}

////////////////////////////////////////////////////////////////////////////////
TFile* RootHistogramsCalib::GetFile(std::string name)  {

  if(pRootFile!=nullptr)
    return pRootFile;
  else{
    std::cout << "Error: Requested file for detector " << name << " does not exist" << std::endl;
    exit(1);
  }
  return 0;
}


