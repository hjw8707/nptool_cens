/*****************************************************************************
 * Copyright (C) 2009-2017   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: morfouace@ganil.fr    *
 *                                                                           *
 * Creation Date  : September 2017                                           *
 * Last update    : November 2024                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Actar Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "TActarPhysics.h"
#include "TActarBeam.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPLNPSystemOfUnits.h
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"

//   ROOT
#include "TChain.h"
#include "TH2F.h"
#include "TCanvas.h"

ClassImp(TActarPhysics)


  ///////////////////////////////////////////////////////////////////////////
TActarPhysics::TActarPhysics()
  : m_EventData(new TActarData),
  m_EventReduced(new MEventReduced),
  m_PreTreatedData(new TActarData),
  m_EventPhysics(this),

  m_Spectra(0),
  fRecoRansac(1),
  fRecoCluster(0),
  fRecoVisu(0),
  fHitThreshold(20),
  fQ_Threshold(0),
  fT_Threshold(0),
  fXBeamMin(0),
  fXBeamMax(128),
  fYBeamMin(60),
  fYBeamMax(67),
  fXTriggerMin1(0),
  fXTriggerMax1(0),
  fYTriggerMin1(0),
  fYTriggerMax1(0),
  fXTriggerMin2(0),
  fXTriggerMax2(0),
  fYTriggerMin2(0),
  fYTriggerMax2(0),
  fAngleBeamMax(5),
  fVertexDistanceMax(5),
  fNumberOfPadsX(128),
  fNumberOfPadsY(128),
  fPadSizeX(2),
  fPadSizeY(2),
  fTimeSampling(80),
  fDriftVelocity(40),
  fPressure(100),
  fGas("iC4H10"),
  fPercentRange(0.95),
  fShortTrackMaxRange(40),
  m_NumberOfPadSilicon(20),
  m_NumberOfDetectors(0),
  fIsGoodEvent(true),
  fScattering(false),
  fIsSimulation(false),
  fPixelStatusFilePath(""),
  fPixelStatusLoaded(false)
  {
    m_randgen = new TRandom3(time(0));
    //InitSpectra();
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TActarPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::AddDetector(double R, double Theta, double Phi, string shape){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos,shape);
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}
///////////////////////////////////////////////////////////////////////////
void TActarPhysics::BuildScatteringPhysicalEvent(){

  if(fRecoRansac && PadX.size()>fHitThreshold){
    m_Track =  m_Ransac->SimpleRansac(PadX, PadY, PadZ, PadCharge);
  }
  else {
    fIsGoodEvent = false;
    return ;
  }

  TrackMult = m_Track.size();
  int voxels_in_all_clusters = 0;


  for(auto & track : m_Track){
    voxels_in_all_clusters+=track.size();
    track.FastLinearFit();
    TVector3 d = track.GetDir();
    double Theta_X = d.Angle(TVector3(1, 0, 0))*180/M_PI;
    ThetaX.push_back(Theta_X);
    if(fabs(Theta_X) < fAngleBeamMax || fabs(180-Theta_X) < fAngleBeamMax){
      BeamTrackMult++;
      if(track.GetDir().X() < 0) track.Inverse();
      m_BeamTrack.push_back(track);
    }
    else{
      ScatteredTrackMult++;
      m_ScatteredTrack.push_back(track);
    }
  }
  if(m_ScatteredTrack.size() == 0){
    fIsGoodEvent = false;
    return ;
  }

  for(auto & bTrack : m_BeamTrack){
    for(auto & sTrack : m_ScatteredTrack){
      TVector3 vertex;
      TVector3 delta;
      double minDist = NPL::MinimumDistanceTwoLines(bTrack.GetP0(), bTrack.GetP0()+bTrack.GetDir(), sTrack.GetP0(), sTrack.GetP0()+sTrack.GetDir(), vertex, delta);
      if(minDist < fVertexDistanceMax){
        if(vertex.X() < 1 || vertex.X() > 255 || vertex.Y() < 110 || vertex.Y() > 146) { //Reject if vertex is out of box and beam zone
          continue;
        }
        sTrack.SetShortTrackMaxLength(fShortTrackMaxRange);
        sTrack.CalculateTrackLength(vertex, fPercentRange); //Do it first cause it also puts the track in the good direction
        if(sTrack.GetFurthestPoint().X() < 3 || sTrack.GetFurthestPoint().X() > 253 || sTrack.GetFurthestPoint().Y() < 3 || sTrack.GetFurthestPoint().Y() > 253){ //Reject if scattered track go out of box
          continue;
        }
        double Phi = TMath::ATan2(sTrack.GetDir().Z(),sTrack.GetDir().Y())>0?TMath::ATan2(sTrack.GetDir().Z(),sTrack.GetDir().Y())*180./M_PI:TMath::ATan2(sTrack.GetDir().Z(),sTrack.GetDir().Y())*180/M_PI+360;
        if(fabs(Phi-90) < 10 || fabs(Phi-270) < 10){ //Reject if track is vertical
          continue;
        }

        TActarScattering aScatteringEvent;
        aScatteringEvent.Vertex = vertex;
        aScatteringEvent.Phi = Phi;
        aScatteringEvent.Range = sTrack.GetTrackLength();
        aScatteringEvent.EndOfTrack = sTrack.GetFurthestPoint();
        aScatteringEvent.Charge = sTrack.GetCharge();
        aScatteringEvent.Theta = bTrack.GetDir().Angle(sTrack.GetDir())*180/M_PI;
        aScatteringEvent.Voxels = sTrack.size();
        std::set<std::pair<int, int>> accountedVoxels; // To keep track of accounted (X, Y) pairs
        for(auto & voxel : sTrack.GetIndex()){
          int x = PadX[voxel];
          int y = PadY[voxel];
          if(IsTriggerZone(x, y) && accountedVoxels.insert({x, y}).second){ // Only count if not already accounted for
            aScatteringEvent.VoxelsInTrigger++;
          }
        }
        TActarBeam aBeam;
        TVector3 beam_entrance(0, bTrack.GetP0().Y()-bTrack.GetDir().Y()/bTrack.GetDir().X()*bTrack.GetP0().X(), bTrack.GetP0().Z()-bTrack.GetDir().Z()/bTrack.GetDir().X()*bTrack.GetP0().X());
        TVector3 beam_exit(256, bTrack.GetP0().Y()-bTrack.GetDir().Y()/bTrack.GetDir().X()*(256-bTrack.GetP0().X()), bTrack.GetP0().Z()-bTrack.GetDir().Z()/bTrack.GetDir().X()*(256-bTrack.GetP0().X()));
        aBeam.SetPoints(beam_entrance, beam_exit);
        aBeam.SetVoxels(bTrack.size());
        aBeam.SetIsScatteredBeam(true);
        Beam.push_back(aBeam);

        VertexMult++;
        //We shorten the track to measure an angle without the deviation
        sTrack.FastLinearFitShort();
        sTrack.CheckDirection(vertex);
        double Phi_Short = TMath::ATan2(sTrack.GetDir().Z(),sTrack.GetDir().Y())>0?TMath::ATan2(sTrack.GetDir().Z(),sTrack.GetDir().Y())*180./M_PI:TMath::ATan2(sTrack.GetDir().Z(),sTrack.GetDir().Y())*180/M_PI+360;

        aScatteringEvent.PhiShort = Phi_Short;
        aScatteringEvent.ThetaShort = bTrack.GetDir().Angle(sTrack.GetDir())*180/M_PI;
        Scattering.push_back(aScatteringEvent);

      }
    }
  }
  if(VertexMult < 1){
    fIsGoodEvent = false;
    return;
  }

}



///////////////////////////////////////////////////////////////////////////
void TActarPhysics::BuildPhysicalEvent() {
  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventReduced = &(**r_ReaderEventData);
  }
  PreTreat();

  //If scattering mode is activated, then we use the dedicated method to process the data
  if(fScattering) {
    BuildScatteringPhysicalEvent();
    return ;
  }
  
  if(fRecoRansac){
    m_Track = m_Ransac->SimpleRansac(PadX, PadY, PadZ, PadCharge);
  }


  //////  NEED TO BE CORRECTED CAUSE TYPE of PADX, PADY changed //////

  // else if(fRecoCluster){
  //   if(PadX.size()>fHitThreshold){
  //     m_Cluster->Init(PadX, PadY, PadZ, PadCharge);
  //     m_Cluster->FindTracks();
  //   }

  //   if(BeamPadX.size()>fHitThreshold){
  //     m_Cluster->Init(BeamPadX, BeamPadY, BeamPadZ, BeamPadCharge);
  //     m_Cluster->FindTracks();
  //   }

  //   m_Track = m_Cluster->GetTracks();

  //   m_Cluster->Clear();
  // }

  TrackMult = m_Track.size();  
  
  //cout << "End of Physical event building" << endl;
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  //ClearPreTreatedData();

  //CleanPads();

  // instantiate CalibrationManager

  //cout << "Start pre-treating" << endl;

  static CalibrationManager* Cal = CalibrationManager::getInstance();


  unsigned int mysize = m_EventReduced->CoboAsad.size();
  for (unsigned int it = 0; it < mysize ; ++it) {
    int co=m_EventReduced->CoboAsad[it].globalchannelid>>11;
    int as=(m_EventReduced->CoboAsad[it].globalchannelid - (co<<11))>>9;
    int ag=(m_EventReduced->CoboAsad[it].globalchannelid - (co<<11)-(as<<9))>>7;
    int ch=m_EventReduced->CoboAsad[it].globalchannelid - (co<<11)-(as<<9)-(ag<<7);
    int where=co*NumberOfASAD*NumberOfAGET*NumberOfChannel + as*NumberOfAGET*NumberOfChannel + ag*NumberOfChannel + ch;
    if(co<16){
      unsigned int vector_size = m_EventReduced->CoboAsad[it].peakheight.size();

      if(fIsSimulation){
        m_EventData->Clear();

        for(unsigned int hh=0; hh<vector_size; hh++){
          int xx = TABLE[4][where], yy = TABLE[5][where], zz = static_cast<int>(m_EventReduced->CoboAsad[it].peaktime[hh]/fTimeSampling+128);
          if(fPixelStatusLoaded && !PixelStatusTable[xx][yy]) continue; //remove the "bad" pixels given by the table
          m_EventData->AddInRebinningMap(1000000*zz+1000*yy+xx, m_EventReduced->CoboAsad[it].peakheight[hh]); //In order to have 512 bin of 80ns in Z
        }
        m_EventData->RebinData();
        for(unsigned int i = 0; i < m_EventData->GetPadMult(); i++){
          int xx = m_EventData->GetPadX(i), yy = m_EventData->GetPadY(i), zz = m_EventData->GetPadZ(i);
          if(fPixelStatusLoaded && !PixelStatusTable[xx][yy]) continue; //remove the "bad" pixels given by the table

          static string name_pixel;
          name_pixel = "Actar/X" ;
          name_pixel+= NPL::itoa(xx);
          name_pixel+= "_Y" ;
          name_pixel+= NPL::itoa(yy);
          name_pixel+= "_Q";
          double CalCharge = CalibrationManager::getInstance()->ApplyCalibration(name_pixel, m_EventData->GetPadCharge(i));
          if(CalCharge>fQ_Threshold){
            double new_zz = (256+zz-m_EventData->GetZTrigger()); //To shift the z according to the trigger time (as in experimental data)
            PadCharge.push_back(CalCharge);
            PadX.push_back((double)xx*fPadSizeX);
            PadY.push_back((double)yy*fPadSizeY);
            PadZ.push_back((double)new_zz*fTimeSampling*fDriftVelocity);
          }
        }
      }
      else{
        for(unsigned int hh=0; hh<vector_size; hh++){
          static string name_pixel;
          name_pixel = "Actar/X" ;
          name_pixel+= NPL::itoa(TABLE[4][where]);
          name_pixel+= "_Y" ;
          name_pixel+= NPL::itoa(TABLE[5][where]);
          name_pixel+= "_Q";
          double CalCharge = CalibrationManager::getInstance()->ApplyCalibration(name_pixel, m_EventReduced->CoboAsad[it].peakheight[hh]);
          if(CalCharge>fQ_Threshold && m_EventReduced->CoboAsad[it].peaktime[hh]>fT_Threshold){
            PadCharge.push_back(CalCharge);
            PadX.push_back((double)TABLE[4][where]*fPadSizeX);
            PadY.push_back((double)TABLE[5][where]*fPadSizeY);
            PadZ.push_back((double)(m_EventReduced->CoboAsad[it].peaktime[hh]*fTimeSampling*fDriftVelocity)); 
          }
        }
      }
    }
    else if(co==31){
      unsigned int vector_size = m_EventReduced->CoboAsad[it].peakheight.size();
      for(unsigned int hit=0;hit<vector_size;hit++){
        if(fInputTreeName=="ACTAR_TTree"){
          int vxi_parameter = m_EventReduced->CoboAsad[it].peaktime[hit];

          if(Si_map[vxi_parameter]<21 && Si_map[vxi_parameter]>0){
            double RawEnergy = m_EventReduced->CoboAsad[it].peakheight[hit];

            static string name;
            name = "ActarSi/D" ;
            name+= NPL::itoa( Si_map[vxi_parameter] - 1) ;
            name+= "_E" ;
            double CalEnergy = CalibrationManager::getInstance()->ApplyCalibration( name, RawEnergy );

            Si_Number.push_back(Si_map[vxi_parameter]);
            Si_E.push_back(CalEnergy);
          }
        }
        else{
          Si_Number.push_back(m_EventReduced->CoboAsad[it].peaktime[hit]);
          Si_E.push_back(m_EventReduced->CoboAsad[it].peakheight[hit]);
        }
      }
    }
  }
  //cout << "End of pre-treatment" << endl;
}

///////////////////////////////////////////////////////////////////////////
bool TActarPhysics::IsBeamZone(int X, int Y)
{
  bool isBeam=false;
  if( (X>=fXBeamMin && X<=fXBeamMax) && (Y>=fYBeamMin && Y<=fYBeamMax) ){
    isBeam=true;
  }

  return isBeam;
}

bool TActarPhysics::IsTriggerZone(int X, int Y){
  bool isTrigger=false;
  if( (X>=fXTriggerMin1 && X<=fXTriggerMax1) && (Y>=fYTriggerMin1 && Y<=fYTriggerMax1) ){
    isTrigger=true;
  }
  else if( (X>=fXTriggerMin2 && X<=fXTriggerMax2) && (Y>=fYTriggerMin2 && Y<=fYTriggerMax2) ){
    isTrigger=true;
  }
  return isTrigger;
}
///////////////////////////////////////////////////////////////////////////
bool TActarPhysics::GoodHit(int iX, int iY)
{
  bool bHit = true;
  if(Hit[iX][iY]<2){
    if(Hit[iX+1][iY]>0 || Hit[iX+1][iY-1]>0 || Hit[iX+1][iY+1]>0){
      if(Hit[iX+2][iY]>0 || Hit[iX+2][iY-1]>0 || Hit[iX+2][iY+1]>0){
        bHit = true;
      }
    }
  }

  return bHit;
}


///////////////////////////////////////////////////////////////////////////
void TActarPhysics::CleanPads()
{
  unsigned int mysize = m_EventReduced->CoboAsad.size();
  for(unsigned int it=0; it<mysize; it++){
    int co=m_EventReduced->CoboAsad[it].globalchannelid>>11;
    int as=(m_EventReduced->CoboAsad[it].globalchannelid - (co<<11))>>9;
    int ag=(m_EventReduced->CoboAsad[it].globalchannelid - (co<<11)-(as<<9))>>7;
    int ch=m_EventReduced->CoboAsad[it].globalchannelid - (co<<11)-(as<<9)-(ag<<7);
    int where=co*NumberOfASAD*NumberOfAGET*NumberOfChannel + as*NumberOfAGET*NumberOfChannel + ag*NumberOfChannel + ch;

    if(co<16){
      unsigned int vector_size = m_EventReduced->CoboAsad[it].peakheight.size();
      for(unsigned int hh=0; hh < vector_size; hh++){
        if(m_EventReduced->CoboAsad[it].peakheight[hh]>fQ_Threshold && m_EventReduced->CoboAsad[it].peaktime[hh]>fT_Threshold){
          Hit[TABLE[4][where]][TABLE[5][where]] += 1;
          /*
             m_PreTreatedData->SetCharge(m_EventReduced->CoboAsad[it].peakheight[hh]);
             m_PreTreatedData->SetPadX(TABLE[4][where]);
             m_PreTreatedData->SetPadY(TABLE[5][where]);
             m_PreTreatedData->SetPadZ(m_EventReduced->CoboAsad[it].peaktime[hh]);*/
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // ACTAR config file //
  string FileName = "./configs/ConfigActar.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigActar.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << "/// Loading user parameter for Analysis from ConfigActar.dat ///" << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigActar.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigActar";
    if (LineBuffer.compare(0, name.length(), name) == 0){
      ReadingStatus = true;
      cout << "**** ConfigActar found" << endl;
    }

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo.compare(0,11,"RecoRansac=") == 0) {
        AnalysisConfigFile >> DataBuffer;
        fRecoRansac = atoi(DataBuffer.c_str());
        cout << "/// Reco using Ransac= " << " " << fRecoRansac << " ///" << endl;
      }

      else if (whatToDo.compare(0,12,"RecoCluster=") == 0) {
        AnalysisConfigFile >> DataBuffer;
        fRecoCluster = atoi(DataBuffer.c_str());
        cout << "/// Reco using Cluster= " << " " << fRecoCluster << " ///" << endl;
      }

      else if (whatToDo.compare(0,9,"RecoVisu=") == 0) {
        AnalysisConfigFile >> DataBuffer;
        fRecoVisu = atoi(DataBuffer.c_str());
        cout << "/// Visualisation= " << " " << fRecoVisu << " ///" << endl;
      }

      else if (whatToDo.compare(0,14,"HIT_THRESHOLD=") == 0) {
        AnalysisConfigFile >> DataBuffer;
        fHitThreshold = atoi(DataBuffer.c_str());
        cout << "/// Hit Threshold= " << " " << fHitThreshold << " ///" << endl;
      }

      else if (whatToDo.compare(0,12,"Q_THRESHOLD=") == 0) {
        AnalysisConfigFile >> DataBuffer;
        fQ_Threshold = atof(DataBuffer.c_str());
        cout << "/// Q Threshold= " << " " << fQ_Threshold << " ///" << endl;
      }

      else if (whatToDo.compare(0,12,"T_THRESHOLD=") == 0) {
        AnalysisConfigFile >> DataBuffer;
        fT_Threshold = atof(DataBuffer.c_str());
        cout << "/// T Threshold= " << " " << fT_Threshold << " ///" << endl;
      }
      else if(whatToDo.compare(0,14,"NumberOfPadsX=")==0){
        AnalysisConfigFile >> DataBuffer;
        fNumberOfPadsX = atoi(DataBuffer.c_str());
        //check_padsX=true;
        cout << "/// Number Of Pads X= " << fNumberOfPadsX << " ///" << endl;
      }

      else if(whatToDo.compare(0,14,"NumberOfPadsY=")==0){
        AnalysisConfigFile >> DataBuffer;
        fNumberOfPadsY = atoi(DataBuffer.c_str());
        //check_padsY=true;
        cout << "/// Number Of Pads Y= " << fNumberOfPadsY << " ///" << endl;
      }

      else if(whatToDo.compare(0,9,"PadSizeX=")==0){
        AnalysisConfigFile >> DataBuffer;
        fPadSizeX = atof(DataBuffer.c_str());
        //check_sizeX=true;
        cout << "/// Pad Size X= " << fPadSizeX << " ///" << endl;
      }

      else if(whatToDo.compare(0,9,"PadSizeY=")==0){
        AnalysisConfigFile >> DataBuffer;
        fPadSizeY = atof(DataBuffer.c_str());
        //check_sizeY=true;
        cout << "/// Pad Size Y= " << fPadSizeY << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"TimeSampling=")==0){
        AnalysisConfigFile >> DataBuffer;
        fTimeSampling = atof(DataBuffer.c_str());
        //check_sizeY=true;
        cout << "/// Time sampling= " << fTimeSampling << " ///" << endl;
      }

      else if(whatToDo.compare(0,9,"XBeamMin=")==0){
        AnalysisConfigFile >> DataBuffer;
        fXBeamMin = atof(DataBuffer.c_str());
        cout << "/// X Beam Min= " << fXBeamMin << " ///" << endl;
      }

      else if(whatToDo.compare(0,9,"XBeamMax=")==0){
        AnalysisConfigFile >> DataBuffer;
        fXBeamMax = atof(DataBuffer.c_str());
        cout << "/// X Beam Max= " << fXBeamMax << " ///" << endl;
      }

      else if(whatToDo.compare(0,9,"YBeamMin=")==0){
        AnalysisConfigFile >> DataBuffer;
        fYBeamMin = atof(DataBuffer.c_str());
        cout << "/// Y Beam Min= " << fYBeamMin << " ///" << endl;
      }

      else if(whatToDo.compare(0,9,"YBeamMax=")==0){
        AnalysisConfigFile >> DataBuffer;
        fYBeamMax = atof(DataBuffer.c_str());
        cout << "/// Y Beam Max= " << fYBeamMax << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"XTriggerMin1=")==0){
        AnalysisConfigFile >> DataBuffer;
        fXTriggerMin1 = atof(DataBuffer.c_str());
        cout << "/// X Trigger Min= " << fXTriggerMin1 << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"XTriggerMax1=")==0){
        AnalysisConfigFile >> DataBuffer;
        fXTriggerMax1 = atof(DataBuffer.c_str());
        cout << "/// X Trigger Max= " << fXTriggerMax1 << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"YTriggerMin1=")==0){
        AnalysisConfigFile >> DataBuffer;
        fYTriggerMin1 = atof(DataBuffer.c_str());
        cout << "/// Y Trigger Min= " << fYTriggerMin1 << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"YTriggerMax1=")==0){
        AnalysisConfigFile >> DataBuffer;
        fYTriggerMax1 = atof(DataBuffer.c_str());
        cout << "/// Y Trigger Max= " << fYTriggerMax1 << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"XTriggerMin2=")==0){
        AnalysisConfigFile >> DataBuffer;
        fXTriggerMin2 = atof(DataBuffer.c_str());
        cout << "/// X Trigger Min= " << fXTriggerMin2 << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"XTriggerMax2=")==0){
        AnalysisConfigFile >> DataBuffer;
        fXTriggerMax2 = atof(DataBuffer.c_str());
        cout << "/// X Trigger Max= " << fXTriggerMax2 << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"YTriggerMin2=")==0){
        AnalysisConfigFile >> DataBuffer;
        fYTriggerMin2 = atof(DataBuffer.c_str());
        cout << "/// Y Trigger Min= " << fYTriggerMin2 << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"YTriggerMax2=")==0){
        AnalysisConfigFile >> DataBuffer;
        fYTriggerMax2 = atof(DataBuffer.c_str());
        cout << "/// Y Trigger Max= " << fYTriggerMax2 << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"AngleBeamMax=")==0){
        AnalysisConfigFile >> DataBuffer;
        fAngleBeamMax = atof(DataBuffer.c_str());
        cout << "/// Angle Beam Max= " << fAngleBeamMax << " ///" << endl;
      }

      else if(whatToDo.compare(0,18,"VertexDistanceMax=")==0){
        AnalysisConfigFile >> DataBuffer;
        fVertexDistanceMax = atof(DataBuffer.c_str());
        cout << "/// Vertex Distance Max= " << fVertexDistanceMax << " ///" << endl;
      }

      else if(whatToDo.compare(0,13,"PercentRange=")==0){
        AnalysisConfigFile >> DataBuffer;
        fPercentRange = atof(DataBuffer.c_str());
        cout << "/// Percentage Range= " << fPercentRange << " ///" << endl;
      }

      else if(whatToDo.compare(0,20,"ShortTrackMaxLength=")==0){
        AnalysisConfigFile >> DataBuffer;
        fShortTrackMaxRange = atof(DataBuffer.c_str());
        cout << "/// Short Track Max Range= " << fShortTrackMaxRange << " ///" << endl;
      }

      else if(whatToDo.compare(0,19,"NumberOfPadSilicon=")==0){
        AnalysisConfigFile >> DataBuffer;
        m_NumberOfPadSilicon = atof(DataBuffer.c_str());
        //check_pressure=true;
        cout << "/// Number of Silicons= " << m_NumberOfPadSilicon << " ///" << endl;
      }

      else if(whatToDo.compare(0,9,"Pressure=")==0){
        AnalysisConfigFile >> DataBuffer;
        fPressure = atof(DataBuffer.c_str());
        //check_pressure=true;
        cout << "/// Pressure= " << fPressure << " ///" << endl;
      }

      else if(whatToDo.compare(0,14,"DriftVelocity=")==0){
        AnalysisConfigFile >> DataBuffer;
        fDriftVelocity = atof(DataBuffer.c_str());
        //check_driftvelocity=true;
        cout << "/// Drift Velocity= " << fDriftVelocity << " ///" << endl;
      }

      else if(whatToDo.compare(0,4,"Gas=")==0){
        AnalysisConfigFile >> DataBuffer;
        fGas = DataBuffer.c_str();
        //check_gas=true;
        cout << "/// Gas Type= " << fGas << " ///" << endl;
      }

      else if(whatToDo.compare(0,11,"Scattering=")==0){
        AnalysisConfigFile >> DataBuffer;
        if(atoi(DataBuffer.c_str()) == 1) {
          fScattering = true;
          cout << "/// Scattering mode activated ///" << endl;
        }
      }

      else if(whatToDo.compare(0,11,"Simulation=")==0){
        AnalysisConfigFile >> DataBuffer;
        if(atoi(DataBuffer.c_str()) == 1) {
          fIsSimulation = true;
          cout << "/// Assuming analysis of simulated data ///" << endl;
        }
      }

      else if(whatToDo.compare(0,12,"PixelStatus=")==0){
        AnalysisConfigFile >> DataBuffer;
          fPixelStatusFilePath = DataBuffer.c_str();
          cout << "/// Pixel Status file path : " << fPixelStatusFilePath << " ///" << endl;
          ifstream pixel_status(fPixelStatusFilePath);
          if(!pixel_status.is_open()) cerr << "Error opening the pixel status file ..." << endl;
          else{
            int xx, yy, vv;
            while(pixel_status >> xx >> yy >> vv){
              PixelStatusTable[xx][yy] = vv>0;
            }
            cout << "Pixel status loaded ! " << endl;
            fPixelStatusLoaded = true;
            pixel_status.close();
          }
      }

    
      else {
        ReadingStatus = false;
      }
    }
      // VXI ACTION FILE //
  string VXI_FileName = "./configs/ACTION_Si_config.dat";
  ifstream VXIConfigFile;
  VXIConfigFile.open(VXI_FileName.c_str());
  if(!VXIConfigFile.is_open()){
    cout << "No VXI ACTION FILE ./configs/ACTION_Si_config.dat found!" << endl;
    //return; //I commented that because I don't use VXI file for my analysis
  }
  else{
    cout << "/// Using VXI ACTION FILE: " << VXI_FileName << " ///" << endl;
    string token;
    int vxi_param, si_nbr;
    for(int i=0; i<20; i++){
      VXIConfigFile >> token >> vxi_param >> si_nbr;
      Si_map[vxi_param] = si_nbr+1;
    }
  }
  VXIConfigFile.close();

  // Lookup table //
  bool ReadingLookupTable = false;
  string LT_FileName; 
  if(IsSimulation()) {LT_FileName = "./configs/LT_simu.dat";}
  else {LT_FileName = "./configs/LT.dat";}
  ifstream LTConfigFile;
  LTConfigFile.open(LT_FileName.c_str());
  if(!LTConfigFile.is_open()){
    cout << "No Lookup Table in ./configs/LT.dat found!" << endl;
    return;
  }
  else{
    cout << "/// Using LookupTable from: " << LT_FileName << " ///" << endl;
    for(int i=0;i<NumberOfCobo*NumberOfASAD*NumberOfAGET*NumberOfChannel;i++){
      LTConfigFile >> TABLE[0][i] >> TABLE[1][i] >> TABLE[2][i] >> TABLE[3][i] >> TABLE[4][i] >> TABLE[5][i];
      //cout << i << "\t" << TABLE[0][i] << "\t" << TABLE[1][i] << "\t" << TABLE[2][i] << "\t" << TABLE[3][i] << "\t" << TABLE[4][i] << "\t" << TABLE[5][i] << endl;
    }
    ReadingLookupTable = true;
  }
  LTConfigFile.close();
  }

  if(fRecoRansac){
    m_Ransac = new NPL::RansacACTAR();
    SetRansacParameter("./configs/RansacConfig_data.dat");
    TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
    asciiConfig->AppendLine("%%% ConfigRANSACActar.dat %%%");
    asciiConfig->Append("./configs/RansacConfig_data.dat");
    asciiConfig->AppendLine("");
  }
  else if(fRecoCluster){
    m_Cluster = new NPL::Cluster(fNumberOfPadsX,fNumberOfPadsY,fRecoVisu);
  }
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::SetRansacParameter(string filename){
  if(fRecoRansac){
    m_Ransac->ReadParameterValue(filename);
  }
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::SetClusterParameter(string filename){
  if(fRecoCluster){
    m_Cluster->ReadParameterValue(filename);
  }
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::Clear() {
  BeamPadX.clear();
  BeamPadY.clear();
  BeamPadZ.clear();
  BeamPadCharge.clear();

  PadX.clear();
  PadY.clear();
  PadZ.clear();
  PadCharge.clear();

  m_Track.clear();
  m_BeamTrack.clear();
  m_ScatteredTrack.clear();

  Si_Number.clear();
  Si_E.clear();

  ThetaX.clear();

  Scattering.clear();
  Beam.clear();

  TrackMult = 0;
  BeamTrackMult = 0;
  ScatteredTrackMult = 0;
  VertexMult = 0;
  fIsGoodEvent = true;
}



///////////////////////////////////////////////////////////////////////////
void TActarPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Actar");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Shape"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Actar " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Actar " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R,Theta,Phi,Shape);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  ReadAnalysisConfig();
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::InitSpectra() {
  m_Spectra = new TActarSpectra(1);
} 

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::SetTreeReader(TTreeReader* TreeReader) {
  TActarPhysicsReader::r_SetTreeReader(TreeReader);
}

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TActarPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TActarPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TActarPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}



////////////////////////////////////////////////////////////////////////////////
vector<TCanvas*> TActarPhysics::GetCanvas() {
  //if(m_Spectra)
  //return m_Spectra->GetCanvas();
  //else{
  vector<TCanvas*> empty;
  return empty;
  //}
}



///////////////////////////////////////////////////////////////////////////
void TActarPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
// void TActarPhysics::AddParameterToCalibrationManager() {
//   CalibrationManager* Cal = CalibrationManager::getInstance();
//   for (int i = 0; i < m_NumberOfDetectors; ++i) {
//     Cal->AddParameter("Actar", "D"+ NPL::itoa(i+1)+"_CHARGE","Actar_D"+ NPL::itoa(i+1)+"_CHARGE");
//     Cal->AddParameter("Actar", "D"+ NPL::itoa(i+1)+"_TIME","Actar_D"+ NPL::itoa(i+1)+"_TIME");
//   }

//   for(int i=0; i<m_NumberOfPadSilicon; i++){
//     Cal->AddParameter("ActarSi","D"+NPL::itoa(i)+"_E","ActarSi_D"+NPL::itoa(i)+"_E");
//   }
// }

void TActarPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int X = 0; X < fNumberOfPadsX; ++X) {
    for (int Y = 0; Y < fNumberOfPadsX; ++Y) {
      Cal->AddParameter("ACTAR", "X"+ NPL::itoa(X)+"_Y"+NPL::itoa(Y)+"_Q","ACTAR_X"+ NPL::itoa(X)+"_Y"+NPL::itoa(Y)+"_Q");
    }
  }
}



///////////////////////////////////////////////////////////////////////////
/*void TActarPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Actar",  true );
  inputChain->SetBranchAddress("Actar", &m_EventData );
  }*/

///////////////////////////////////////////////////////////////////////////
void TActarPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();

  // Option to use the npreader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
    inputChain->SetBranchStatus("data",  true );
    //inputChain->SetBranchStatus("ACTAR_TTree",  true );
    inputChain->SetBranchAddress("data", &m_EventReduced );
  }
  fInputTreeName = inputChain->GetTree()->GetName();
}



///////////////////////////////////////////////////////////////////////////
void TActarPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
    inputChain->SetBranchAddress("Actar", &m_EventPhysics);
  }
}



///////////////////////////////////////////////////////////////////////////
void TActarPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Actar", "TActarPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TActarPhysics::Construct() { return (NPL::VDetector*) new TActarPhysics();}

NPL::VTreeReader* TActarPhysics::ConstructReader() { return (NPL::VTreeReader*)new TActarPhysicsReader();}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Actar{
  public:
    proxy_Actar(){
      NPL::DetectorFactory::getInstance()->AddToken("Actar","Actar");
      NPL::DetectorFactory::getInstance()->AddDetector("Actar",TActarPhysics::Construct);
      NPL::DetectorFactory::getInstance()->AddDetectorReader("Actar", TActarPhysics::ConstructReader);

    }
};

proxy_Actar p_Actar;
}
