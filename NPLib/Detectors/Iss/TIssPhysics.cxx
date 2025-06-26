/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 * Author: M. Labiche                     address: marc.labiche@stfc.ac.uk   *
 *                                                                           *
 * Creation Date  : Jul 2019                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Iss treated data                                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TIssPhysics.h"
using namespace ISS_LOCAL;

//   STL
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdlib.h>

//   NPL
#include "NPDetectorFactory.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"

using namespace NPUNITS;
//   ROOT
#include "TChain.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TIssPhysics)

///////////////////////////////////////////////////////////////////////////
TIssPhysics::TIssPhysics() {
	EventMultiplicity                  = 0;
	m_EventData                        = new TIssData;
	m_PreTreatedData                   = new TIssData;
	m_EventPhysics                     = this;
	m_Spectra                          = NULL;
	//  m_NumberOfTelescope                = 0;
	m_NumberOfDetector                = 0;
	
	m_MaximumStripMultiplicityAllowed  = 10;
	m_StripEnergyMatchingSigma         = 0.020;
	m_StripEnergyMatchingNumberOfSigma = 3;
	
	// Threshold
	m_StripFront_E_RAW_Threshold = 0 ;
	m_StripFront_E_Threshold = 0 ;
	
	m_StripBack_E_RAW_Threshold = 0 ;
	m_StripBack_E_Threshold = 0 ;
	
	m_Take_E_Front=true;
	m_Take_T_Back=true;
	
	m_NominalField = 0;
	
	
	/* from Must2
	 // Raw Threshold
	 m_Si_X_E_RAW_Threshold = 8200;
	 m_Si_Y_E_RAW_Threshold = 8200;
	 // Calibrated Threshold
	 m_Si_X_E_Threshold = 0;
	 m_Si_Y_E_Threshold = 0;
	 
	 
	 m_Take_E_Y = false;
	 m_Take_T_Y = true;
	 */
	//////////////////
	
}

///////////////////////////////////////////////////////////////////////////
TIssPhysics::~TIssPhysics() {}



///////////////////////////////////////////////////////////////////////////

void TIssPhysics::BuildPhysicalEvent() {
	
	PreTreat();
	
	vector< TVector2 > couple = Match_Front_Back() ;
	EventMultiplicity = couple.size();
	
	
	/*
	 for(unsigned int i = 0 ; i < EventMultiplicity ; ++i){
	 cout << " iX= " << couple[i].X() << "  iY=" << couple[i].Y() << endl;
	 }
	 //CHECK
	 */
	
	unsigned int size = couple.size();
	for(unsigned int i = 0 ; i < size ; ++i){
		
		int N = m_PreTreatedData->GetFront_DetectorNbr(couple[i].X()) ;
		
		int Front = m_PreTreatedData->GetFront_StripNbr(couple[i].X()) ;
		int Back  = m_PreTreatedData->GetBack_StripNbr(couple[i].Y()) ;
		
		double Front_E = m_PreTreatedData->GetFront_Energy( couple[i].X() ) ;
		double Back_E  = m_PreTreatedData->GetBack_Energy( couple[i].Y() ) ;
		
		double Front_T = m_PreTreatedData->GetFront_TimeCFD( couple[i].X() ) ;
		double Back_T  = m_PreTreatedData->GetBack_TimeCFD ( couple[i].Y() ) ;
		
		DetectorNumber.push_back(N);
		StripFront_E.push_back(Front_E);
		StripFront_T.push_back(Front_T) ;
		StripBack_E.push_back(Back_E) ;
		StripBack_T.push_back(Back_T) ;
		Strip_Front.push_back(Front) ;
		Strip_Back.push_back(Back) ;
		
		// Try to obtain Pixel Calibration
		static CalibrationManager* Cal = CalibrationManager::getInstance();
		static string name;
		//Store for calibration purposes
		Strip_Front_RawE.push_back(StripFront_OriginalE[couple[i].X()]);
		Strip_Back_RawE.push_back(StripBack_OriginalE[couple[i].Y()]);
		
		//Proceed for Pixel Calibration
		name = "ISS/D"+ NPL::itoa(N)+"_STRIP_FRONT"+ NPL::itoa(Front)+"_BACK"+ NPL::itoa(Back)+"_E";
		double Pixel_E = Cal->ApplyCalibration(name,StripFront_OriginalE[couple[i].X()] );
		if(Pixel_E != StripFront_OriginalE[couple[i].X()]){
			Strip_E.push_back(Pixel_E);
			name = "ISS/D"+ NPL::itoa(N)+"_STRIP_FRONT"+ NPL::itoa(Front)+"_BACK"+ NPL::itoa(Back)+"_DEADLAYER";
			DeadLayer.push_back(Cal->GetValue(name,0));
		}
		// Fall Back option, take the Strip Calibration
		else if(m_Take_E_Front){
			Strip_E.push_back(Front_E) ;
			name = "ISS/D"+ NPL::itoa(N)+"_STRIP_FRONT"+ NPL::itoa(Front)+"_DEADLAYER";
			DeadLayer.push_back(Cal->GetValue(name,0));
		}
		else{
			Strip_E.push_back(Back_E) ;
			name = "ISS/D"+ NPL::itoa(N)+"_STRIP_BACK"+ NPL::itoa(Back)+"_DEADLAYER";
			DeadLayer.push_back(Cal->GetValue(name,0));
		}
		
		if(m_Take_T_Back)
			Strip_T.push_back(Back_T) ;
		else
			Strip_T.push_back(Front_T) ;
		
		
	}
	
	if(DetectorNumber.size()==1)
		return;
	
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TIssPhysics::BuildSimplePhysicalEvent(){ 
	
	
	BuildPhysicalEvent();
	
}



///////////////////////////////////////////////////////////////////////////
void TIssPhysics::PreTreat() {
	
	
	ClearPreTreatedData();
	
	
	//   Front
	unsigned int sizeFront = m_EventData->GetMultiplicityFront();
	for(unsigned int i = 0 ; i < sizeFront ; i++){
		if( m_EventData->GetFront_Energy(i)>m_StripFront_E_RAW_Threshold && IsValidChannel("Front", m_EventData->GetFront_DetectorNbr(i), m_EventData->GetFront_StripNbr(i)) ){
			double Front_E = fStrip_Front_E(m_EventData , i);
			if( Front_E > m_StripFront_E_Threshold ){
				m_PreTreatedData->SetFront( m_EventData->GetFront_DetectorNbr(i),
										   m_EventData->GetFront_StripNbr(i),
										   Front_E,
										   m_EventData->GetFront_TimeCFD(i),
										   m_EventData->GetFront_TimeLED(i));
				
				StripFront_OriginalE.push_back( m_EventData->GetFront_Energy(i));
			}
		}
	}
	
	//  Back
	unsigned int sizeBack = m_EventData->GetMultiplicityBack() ;
	for(unsigned int i = 0 ; i < sizeBack ; i++){
		if( m_EventData->GetBack_Energy(i)>m_StripBack_E_RAW_Threshold && IsValidChannel("Back", m_EventData->GetBack_DetectorNbr(i), m_EventData->GetBack_StripNbr(i)) ){
			double Back_E = fStrip_Back_E(m_EventData , i);
			if( Back_E > m_StripBack_E_Threshold ){
				m_PreTreatedData->SetBack( m_EventData->GetBack_DetectorNbr(i),
										  m_EventData->GetBack_StripNbr(i),
										  Back_E,
										  m_EventData->GetBack_TimeCFD(i),
										  m_EventData->GetBack_TimeLED(i));
				
				StripBack_OriginalE.push_back( m_EventData->GetBack_Energy(i));
			}
		}
	}
	
	
	
	return;
}


///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
int TIssPhysics :: CheckEvent(){
	return 1 ; // Regular Event
}

///////////////////////////////////////////////////////////////////////////



vector < TVector2 > TIssPhysics :: Match_Front_Back(){
	vector < TVector2 > ArrayOfGoodCouple ;
	// Prevent code from treating very high multiplicity Event
	// Those event are not physical anyway and that improve speed.
	if( m_PreTreatedData->GetMultiplicityFront() > m_MaximumStripMultiplicityAllowed || m_PreTreatedData->GetMultiplicityBack() > m_MaximumStripMultiplicityAllowed )
		return ArrayOfGoodCouple;
	
	ArrayOfGoodCouple.clear();
	unsigned int mysizeF =  m_PreTreatedData->GetMultiplicityFront();
	unsigned int mysizeB =  m_PreTreatedData->GetMultiplicityBack();
	//cout << " mysizeF/B = " << mysizeF << " " << mysizeB << endl; //CHECK
	
	for(unsigned int i = 0 ; i < mysizeF ; i++) {
		for(unsigned int j = 0 ; j < mysizeB ; j++){
			//cout << "Det Number F/B = " << m_PreTreatedData->GetFront_DetectorNbr(i) << " =? " << m_PreTreatedData->GetBack_DetectorNbr(j) << endl;
			//   if same detector check energy
			if ( m_PreTreatedData->GetFront_DetectorNbr(i) == m_PreTreatedData->GetBack_DetectorNbr(j) ){
				//cout << "Energy F= " << m_PreTreatedData->GetFront_Energy(i) << endl;
				//cout << "Energy B= " << m_PreTreatedData->GetBack_Energy(j) << endl;
				//cout << "abs(F-B)= " << abs( (m_PreTreatedData->GetFront_Energy(i)-m_PreTreatedData->GetBack_Energy(j))/2. ) << endl;
				//cout << "Allowed difference = " << m_StripEnergyMatchingNumberOfSigma*m_StripEnergyMatchingSigma << endl;
				
				//   Look if energy match
				if( abs( (m_PreTreatedData->GetFront_Energy(i)-m_PreTreatedData->GetBack_Energy(j))/2. ) < m_StripEnergyMatchingNumberOfSigma*m_StripEnergyMatchingSigma ){
					//cout << " Good Couple!! " << endl;
					ArrayOfGoodCouple.push_back ( TVector2(i,j) ) ;
				}
				
			}
		}
	}
	
	//cout << " Good couples = " << ArrayOfGoodCouple.size()  << "  must.be.less.then  " << m_PreTreatedData->GetMultiplicityFront() << endl;
	
	//  Prevent to treat event with ambiguous matching beetween Front and Back
	if( ArrayOfGoodCouple.size() > m_PreTreatedData->GetMultiplicityFront() )
		ArrayOfGoodCouple.clear() ;
	
	//cout << " final size = " << ArrayOfGoodCouple.size() << endl;
	
	return ArrayOfGoodCouple;
	
	
	
	
	
	/* from MUST2
	 vector<TVector2> TIssPhysics::Match_X_Y() {
	 vector<TVector2> ArrayOfGoodCouple;
	 ArrayOfGoodCouple.clear();
	 
	 m_StripXEMult = m_PreTreatedData->GetMMStripXEMult();
	 m_StripYEMult = m_PreTreatedData->GetMMStripYEMult();
	 
	 double matchSigma  = m_StripEnergyMatchingSigma;
	 double NmatchSigma = m_StripEnergyMatchingNumberOfSigma;
	 
	 // Prevent code from treating very high multiplicity Event
	 // Those event are not physical anyway and that improve speed.
	 if (m_StripXEMult > m_MaximumStripMultiplicityAllowed
	 || m_StripYEMult > m_MaximumStripMultiplicityAllowed) {
	 return ArrayOfGoodCouple;
	 }
	 
	 // Get Detector multiplicity
	 for (unsigned int i = 0; i < m_StripXEMult; i++) {
	 int N = m_PreTreatedData->GetMMStripXEDetectorNbr(i);
	 m_StripXMultDet[N] += 1;
	 }
	 
	 for (unsigned int j = 0; j < m_StripYEMult; j++) {
	 int N = m_PreTreatedData->GetMMStripYEDetectorNbr(j);
	 m_StripYMultDet[N] += 1;
	 }
	 
	 for (unsigned int i = 0; i < m_StripXEMult; i++) {
	 for (unsigned int j = 0; j < m_StripYEMult; j++) {
	 
	 // Declaration of variable for clarity
	 int StripXDetNbr = m_PreTreatedData->GetMMStripXEDetectorNbr(i);
	 int StripYDetNbr = m_PreTreatedData->GetMMStripYEDetectorNbr(j);
	 
	 //   if same detector check energy
	 if (StripXDetNbr == StripYDetNbr) {
	 
	 int DetNbr = StripXDetNbr;
	 
	 // Declaration of variable for clarity
	 double StripXEnergy = m_PreTreatedData->GetMMStripXEEnergy(i);
	 double StripXNbr    = m_PreTreatedData->GetMMStripXEStripNbr(i);
	 
	 double StripYEnergy = m_PreTreatedData->GetMMStripYEEnergy(j);
	 double StripYNbr    = m_PreTreatedData->GetMMStripYEStripNbr(j);
	 
	 //   Look if energy match
	 // FIXME Should be proportional to the energy loss in the DSSDs
	 // if (abs(StripXEnergy - StripYEnergy)
	 // < 0.09 * (std::max(StripXEnergy, StripYEnergy))) {
	 // negligible correction according to Adrien
	 /*
	 if (abs((StripXEnergy - StripYEnergy) / 2.)
	 < NmatchSigma * matchSigma) {
	 
	 // Special Option, if the event is between two CsI
	 // cristal, it is rejected.
	 if (m_Ignore_not_matching_CsI) {
	 bool check_validity = false;
	 for (unsigned int hh = 0; hh < 16; ++hh) {
	 if (Match_Si_CsI(StripXNbr, StripYNbr, hh + 1)) {
	 check_validity = true;
	 }
	 }
	 if (check_validity) {
	 ArrayOfGoodCouple.push_back(TVector2(i, j));
	 }
	 }
	 
	 // Special Option, if the event is between two SiLi pad ,
	 // it is rejected.
	 else if (m_Ignore_not_matching_SiLi) {
	 bool check_validity = false;
	 for (unsigned int hh = 0; hh < 16; ++hh) {
	 if (Match_Si_SiLi(StripXNbr, StripYNbr, hh + 1))
	 check_validity = true;
	 }
	 if (check_validity)
	 ArrayOfGoodCouple.push_back(TVector2(i, j));
	 }
	 // Regular case, keep the event
	 else {
	 ArrayOfGoodCouple.push_back(TVector2(i, j));
	 m_NMatchDet[DetNbr] += 1;
	 }
	 }
	 } // if same detector
	 } // loop on StripY Mult
	 } // loop on StripX Mult
	 
	 unsigned int couple_size = ArrayOfGoodCouple.size();
	 
	 for (unsigned int i = 0; i < couple_size; ++i) {
	 int N = m_PreTreatedData->GetMMStripXEDetectorNbr(ArrayOfGoodCouple[i].X());
	 m_match_type.push_back(CheckEvent(N));
	 }
	 return ArrayOfGoodCouple;
	 */
	
}



////////////////////////////////////////////////////////////////////////////
bool TIssPhysics :: IsValidChannel(const string& DetectorType, const int& telescope , const int& channel){
	
	if(DetectorType == "Front")
		return *(m_FrontChannelStatus[telescope-1].begin()+channel-1);
	
	else if(DetectorType == "Back")
		return *(m_BackChannelStatus[telescope-1].begin()+channel-1);
	
	else return false;
}

///////////////////////////////////////////////////////////////////////////
void TIssPhysics::ReadAnalysisConfig(){
	bool ReadingStatus = false;
	
	// path to file
	string FileName = "./configs/ConfigIss.dat";
	
	// open analysis config file
	ifstream AnalysisConfigFile;
	AnalysisConfigFile.open(FileName.c_str());
	
	if (!AnalysisConfigFile.is_open()) {
		cout << " No ConfigIss.dat found: Default parameter loaded for Analysis " << FileName << endl;
		return;
	}
	cout << " Loading user parameter for Analysis from ConfigIss.dat " << endl;
	
	// Save it in a TAsciiFile
	TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
	asciiConfig->AppendLine("%%% ConfigIss.dat %%%");
	asciiConfig->Append(FileName.c_str());
	asciiConfig->AppendLine("");
	// read analysis config file
	string LineBuffer,DataBuffer,whatToDo;
	while (!AnalysisConfigFile.eof()) {
		// Pick-up next line
		getline(AnalysisConfigFile, LineBuffer);
		
		// search for "header"
		if (LineBuffer.compare(0, 11, "ConfigIss") == 0) ReadingStatus = true;
		
		// loop on tokens and data
		while (ReadingStatus ) {
			
			whatToDo="";
			AnalysisConfigFile >> whatToDo;
			
			// Search for comment symbol (%)
			if (whatToDo.compare(0, 1, "%") == 0) {
				AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
			}
			
			else if (whatToDo=="MAX_STRIP_MULTIPLICITY") {
				AnalysisConfigFile >> DataBuffer;
				m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str() );
				cout << "MAXIMUN STRIP MULTIPLICITY " << m_MaximumStripMultiplicityAllowed << endl;
			}
			
			else if (whatToDo=="STRIP_ENERGY_MATCHING_SIGMA") {
				AnalysisConfigFile >> DataBuffer;
				m_StripEnergyMatchingSigma = atof(DataBuffer.c_str() );
				cout << "STRIP ENERGY MATCHING SIGMA " << m_StripEnergyMatchingSigma <<endl;
			}
			
			else if (whatToDo=="STRIP_ENERGY_MATCHING_NUMBER_OF_SIGMA") {
				AnalysisConfigFile >> DataBuffer;
				m_StripEnergyMatchingNumberOfSigma = atoi(DataBuffer.c_str() );
				cout << "STRIP ENERGY MATCHING NUMBER OF SIGMA " << m_StripEnergyMatchingNumberOfSigma << endl;
			}
			
			else if (whatToDo== "DISABLE_ALL") {
				AnalysisConfigFile >> DataBuffer;
				cout << whatToDo << "  " << DataBuffer << endl;
				int Detector = atoi(DataBuffer.substr(2,2).c_str());
				vector< bool > ChannelStatus;
				ChannelStatus.resize(128,false);
				//ChannelStatus.resize(24,false);
				m_FrontChannelStatus[Detector-1] = ChannelStatus;
				ChannelStatus.resize(11,false);
				//ChannelStatus.resize(22,false);
				//ChannelStatus.resize(48,false);
				m_BackChannelStatus[Detector-1] = ChannelStatus;
				
			}
			
			else if (whatToDo == "DISABLE_CHANNEL") {
				AnalysisConfigFile >> DataBuffer;
				int Detector = atoi(DataBuffer.substr(2,2).c_str());
				int channel = -1;
				if (DataBuffer.find("STRF") != string::npos) {
					channel = atoi(DataBuffer.substr(8).c_str());
					*(m_FrontChannelStatus[Detector-1].begin()+channel-1) = false;
					cout << "DISABLE DETECTOR " << Detector << " STRIP FRONT " << channel << endl;
				}
				
				else if (DataBuffer.find("STRB")!=string::npos) {
					channel = atoi(DataBuffer.substr(8).c_str());
					*(m_BackChannelStatus[Detector-1].begin()+channel-1) = false;
					cout << "DISABLE DETECTOR " << Detector << " STRIP BACK " << channel << endl;
					
				}
				
				
				else cout << "Warning: detector type for Iss unknown!" << endl;
				
			}
			
			else if (whatToDo=="TAKE_E_FRONT") {
				m_Take_E_Front = true;
				cout << whatToDo << endl;
			}
			
			else if (whatToDo=="TAKE_E_BACK") {
				m_Take_E_Front = false;
				cout << whatToDo << endl;
			}
			
			else if (whatToDo=="TAKE_T_FRONT") {
				m_Take_T_Back = false;
				cout << whatToDo << endl;
			}
			
			else if (whatToDo=="TAKE_T_BACK") {
				m_Take_T_Back = true;
				cout << whatToDo << endl;
			}
			
			else if (whatToDo=="STRIP_FRONT_E_RAW_THRESHOLD") {
				AnalysisConfigFile >> DataBuffer;
				m_StripFront_E_RAW_Threshold = atof(DataBuffer.c_str());
				cout << whatToDo << " " << m_StripFront_E_RAW_Threshold << endl;
			}
			
			else if (whatToDo=="STRIP_BACK_E_RAW_THRESHOLD") {
				AnalysisConfigFile >> DataBuffer;
				m_StripBack_E_RAW_Threshold = atof(DataBuffer.c_str());
				cout << whatToDo << " " << m_StripBack_E_RAW_Threshold << endl;
			}
			
			else if (whatToDo=="STRIP_FRONT_E_THRESHOLD") {
				AnalysisConfigFile >> DataBuffer;
				m_StripFront_E_Threshold = atof(DataBuffer.c_str());
				cout << whatToDo << " " << m_StripFront_E_Threshold << endl;
			}
			
			else if (whatToDo=="STRIP_BACK_E_THRESHOLD") {
				AnalysisConfigFile >> DataBuffer;
				m_StripBack_E_Threshold = atof(DataBuffer.c_str());
				cout << whatToDo << " " << m_StripBack_E_Threshold << endl;
			}
			
			
			else {
				ReadingStatus = false;
			}
			
		}
	}
	
}


///////////////////////////////////////////////////////////////////////////
void TIssPhysics::Clear(){
	EventMultiplicity = 0;
	
	//   Provide a Classification of Event
	EventType.clear() ;
	
	// Detector
	DetectorNumber.clear() ;
	
	//   DSSD
	Strip_E.clear() ;
	Strip_T.clear() ;
	StripFront_E.clear() ;
	StripFront_T.clear();
	StripBack_E.clear() ;
	StripBack_T.clear() ;
	Strip_Front.clear() ;
	Strip_Back.clear() ;
	StripFront_OriginalE.clear();
	StripBack_OriginalE.clear();
	DeadLayer.clear();
	Strip_Front_RawE.clear();
	Strip_Back_RawE.clear();
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TIssPhysics::ReadConfiguration(NPL::InputParser parser){
	vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ISSDet");
	if(NPOptionManager::getInstance()->GetVerboseLevel())
		cout << "//// " << blocks.size() << " detectors found " << endl;
	
	// Cartesian Case
	vector<string> cart
	= {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128", "SI"};
	vector<string> cart1
	= {"MField", "X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128", "SI"};
	
	// Spherical Case
	vector<string> sphe = {"R", "THETA", "PHI", "BETA", "SI"};
	
	for (unsigned int i = 0; i < blocks.size(); i++) {
		
		if (blocks[i]->HasTokenList(cart1)) {
			if (NPOptionManager::getInstance()->GetVerboseLevel())
				cout << endl << "////  Iss Telescope " << i + 1 << endl;
			m_NominalField= blocks[i]->GetDouble("MField", "tesla");
			TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
			TVector3 B = blocks[i]->GetTVector3("X128_Y1", "mm");
			TVector3 C = blocks[i]->GetTVector3("X1_Y128", "mm");
			TVector3 D = blocks[i]->GetTVector3("X128_Y128", "mm");
			AddTelescope(A, B, C, D);
		}
		
		else if (blocks[i]->HasTokenList(cart)) {
			if (NPOptionManager::getInstance()->GetVerboseLevel())
				cout << endl << "////  Iss Telescope " << i + 1 << endl;
			TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
			TVector3 B = blocks[i]->GetTVector3("X128_Y1", "mm");
			TVector3 C = blocks[i]->GetTVector3("X1_Y128", "mm");
			TVector3 D = blocks[i]->GetTVector3("X128_Y128", "mm");
			AddTelescope(A, B, C, D);
			
		}
		
		else if (blocks[i]->HasTokenList(sphe)) {
			if (NPOptionManager::getInstance()->GetVerboseLevel())
				cout << endl << "////  Iss Telescope " << i + 1 << endl;
			
			double         Theta = blocks[i]->GetDouble("THETA", "deg");
			double         Phi   = blocks[i]->GetDouble("PHI", "deg");
			double         R     = blocks[i]->GetDouble("R", "mm");
			vector<double> beta  = blocks[i]->GetVectorDouble("BETA", "deg");
			//AddTelescope(Theta, Phi, R, beta[0], beta[1], beta[2]);
			
		}
		
		else {
			cout << "ERROR: Missing token for M2Telescope blocks, check your "
			"input "
			"file"
			<< endl;
			exit(1);
		}
	}
	
	
	/*
	 vector<string> tokenBOX =     {"Z","ThicknessDetector1","ThicknessDetector2","ThicknessDetector3","ThicknessDetector4","ThicknessPAD1","ThicknessPAD2","ThicknessPAD3","ThicknessPAD4"};
	 
	 for(unsigned int i = 0 ; i < blocks.size() ; i++){
	 
	 if(blocks[i]->GetMainValue()=="BOX" && blocks[i]->HasTokenList(tokenBOX)){
	 if(NPOptionManager::getInstance()->GetVerboseLevel())
	 cout << endl << "////  Iss Box " << i+1 <<  endl;
	 double Z = blocks[i]->GetDouble("Z","mm");
	 double Thickness1= blocks[i]->GetDouble("ThicknessDetector1","micrometer");
	 double Thickness2= blocks[i]->GetDouble("ThicknessDetector2","micrometer");
	 double Thickness3= blocks[i]->GetDouble("ThicknessDetector3","micrometer");
	 double Thickness4= blocks[i]->GetDouble("ThicknessDetector4","micrometer");
	 double ThicknessPAD1 = blocks[i]->GetDouble("ThicknessPAD1","micrometer");
	 double ThicknessPAD2 = blocks[i]->GetDouble("ThicknessPAD2","micrometer");
	 double ThicknessPAD3 = blocks[i]->GetDouble("ThicknessPAD3","micrometer");
	 double ThicknessPAD4 = blocks[i]->GetDouble("ThicknessPAD4","micrometer");
	 AddBoxDetector(Z);
	 }
	 
	 else{
	 cout << "Warning: check your input file formatting " << endl;
	 }
	 }
	 */
	
	InitializeStandardParameter();
	ReadAnalysisConfig();
}
///////////////////////////////////////////////////////////////////////////
void TIssPhysics::InitSpectra(){  
	m_Spectra = new TIssSpectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TIssPhysics::FillSpectra(){  
	m_Spectra -> FillRawSpectra(m_EventData);
	m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
	m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TIssPhysics::CheckSpectra(){  
	m_Spectra->CheckSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TIssPhysics::ClearSpectra(){  
	// To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TIssPhysics::GetSpectra() {
	if(m_Spectra)
		return m_Spectra->GetMapHisto();
	else{
		map< string , TH1*> empty;
		return empty;
	}
} 

///////////////////////////////////////////////////////////////////////////
void TIssPhysics::WriteSpectra(){
	m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TIssPhysics::AddParameterToCalibrationManager(){
	CalibrationManager* Cal = CalibrationManager::getInstance();
	
	for(int i = 0 ; i < m_NumberOfDetector ; ++i){
		
		for( int j = 0 ; j < 128 ; ++j){
			// Front Strip Calibration
			Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_E","ISS_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_E")   ;
			Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_DEADLAYER",
							  "ISS_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_DEADLAYER")   ;
			Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_T","ISS_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_T")   ;
			
			// Pixel Calibration
			//for( int k = 0 ; k < 22 ; ++k){
			for( int k = 0 ; k < 11 ; ++k){
				
				Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_BACK"+ NPL::itoa(k+1)+"_E","ISS_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_BACK"+ NPL::itoa(k+1)+"_E")   ;
				Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_BACK"+ NPL::itoa(k+1)+"_DEADLAYER",
								  "ISS_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_BACK"+ NPL::itoa(k+1)+"_DEADLAYER")   ;
				
			}
		}
		
		//for( int j = 0 ; j < 22 ; ++j){
		for( int j = 0 ; j < 11 ; ++j){
			
			// Back strip Calibration
			Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_E","ISS_D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_E")   ;
			Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_DEADLAYER",
							  "ISS_D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_DEADLAYER")   ;
			
			Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_T","ISS_D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_T")   ;
		}
		
		for( int j = 0 ; j < 1 ; ++j){
			// Pad Calibration
			Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_PAD_E","ISS_D"+ NPL::itoa(i+1)+"_PAD_E")   ;
			Cal->AddParameter("ISS", "D"+ NPL::itoa(i+1)+"_PAD_T","ISS_D"+ NPL::itoa(i+1)+"_PAD_T")   ;
		}
	}
	
	return;
	
}

///////////////////////////////////////////////////////////////////////////
void TIssPhysics::InitializeRootInputRaw(){
	TChain* inputChain = RootInput::getInstance()->GetChain()   ;
	inputChain->SetBranchStatus( "ISS" , true );
	// The following line is necessary only for system were the tree is splitted
	// (older root version). The found argument silenced the Branches not found
	// warning for non splitted tree.
	if(inputChain->FindBranch("fIss_*"))
		inputChain->SetBranchStatus( "fIss_*",true);
	inputChain->SetBranchAddress( "ISS" , &m_EventData );
	
}

///////////////////////////////////////////////////////////////////////////
void TIssPhysics::InitializeRootInputPhysics(){
	TChain* inputChain = RootInput::getInstance()->GetChain();
	inputChain->SetBranchStatus( "EventMultiplicity" , true );
	inputChain->SetBranchStatus( "EventType" , true );
	inputChain->SetBranchStatus( "DetectorNumber" , true );
	inputChain->SetBranchStatus( "Strip_E" , true );
	inputChain->SetBranchStatus( "Strip_T" , true );
	inputChain->SetBranchStatus( "StripFront_E" , true );
	inputChain->SetBranchStatus( "StripFront_T" , true );
	inputChain->SetBranchStatus( "StripBack_E" , true );
	inputChain->SetBranchStatus( "StripBack_T" , true );
	inputChain->SetBranchStatus( "Strip_Front" , true );
	inputChain->SetBranchStatus( "Strip_Back" , true );
	inputChain->SetBranchAddress( "ISS" , &m_EventPhysics )      ;
}

///////////////////////////////////////////////////////////////////////////
void TIssPhysics::InitializeRootOutput(){
	TTree* outputTree = RootOutput::getInstance()->GetTree();
	outputTree->Branch( "ISS" , "TIssPhysics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
/////   Specific to IssArray   ////

// telescope
void TIssPhysics::AddTelescope(TVector3 C_X1_Y1, TVector3 C_X128_Y1,
							   TVector3 C_X1_Y128, TVector3 C_X128_Y128) {
	// To avoid warning
	C_X128_Y128 *= 1;
	
	//m_NumberOfTelescope++;
	m_NumberOfDetector++;
	
	//   Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that
	//   Y strip are allong X axis)
	TVector3 U      = C_X128_Y1 - C_X1_Y1;
	//double   Ushift = (U.Mag() - 98) / 2.; // square Si
	//double   Ushift = (U.Mag() - 22) / 2.;
	U               = U.Unit();
	//   Vector V on Telescope Face (parallele to X Strip)
	TVector3 V      = C_X1_Y128 - C_X1_Y1;
	//  double   Vshift = (V.Mag() - 98) / 2.; // square Si
	//double   Vshift = (V.Mag() - 125) / 2.;
	V               = V.Unit();
	
	//   Position Vector of Strip Center
	TVector3 StripCenter = TVector3(0, 0, 0);
	//   Position Vector of X=1 Y=1 Strip
	TVector3 Strip_1_1;
	
	//   Geometry Parameter
	//double Face          = 98; // mm
	//double NumberOfStrip = 128;
	//double StripPitch    = Face / NumberOfStrip; // mm
	
	double FaceL          = 22; // mm
	double NumberOfLStrip = 11;
	double StripLPitch    = FaceL / NumberOfLStrip; // mm
	
	double FaceT          = 125; // mm
	double NumberOfTStrip = 128;
	double StripTPitch    = FaceT / NumberOfTStrip; // mm
	
	
	//   Buffer object to fill Position Array
	vector<double> lineX;
	vector<double> lineY;
	vector<double> lineZ;
	
	vector<vector<double>> OneTelescopeStripPositionX;
	vector<vector<double>> OneTelescopeStripPositionY;
	vector<vector<double>> OneTelescopeStripPositionZ;
	
	//   Moving StripCenter to 1.1 corner:
	//Strip_1_1 = C_X1_Y1 + (U + V) * (StripPitch / 2.);
	//Strip_1_1 += U * Ushift + V * Vshift;
	
	Strip_1_1 = C_X1_Y1 + (U*(StripLPitch/2.)) + (V*(StripTPitch/2.));  // or: C_X1_Y1 + (U * (StripTPitch/2.) + (V* (StripLPitch/2.) !!!  To be checked !!!
	
	cout << "Strip_1_1 x=" << Strip_1_1[0] << endl;
	cout << "Strip_1_1 y=" << Strip_1_1[1] << endl;
	cout << "Strip_1_1 z=" << Strip_1_1[2] << endl;
	
	//Strip_1_1 += U * Ushift + V * Vshift; //
	
	//cout << "Shifted Strip_1_1 x=" << Strip_1_1[0] << endl;
	//cout << "Shifted Strip_1_1 y=" << Strip_1_1[1] << endl;
	//cout << "Shifted Strip_1_1 z=" << Strip_1_1[2] << endl;
	
	
	//for (int i = 0; i < 128; ++i) {
	for (int i = 0; i < NumberOfTStrip; ++i) {  // front strip are the transversal one
		lineX.clear();
		lineY.clear();
		lineZ.clear();
		
		//for (int j = 0; j < 128; ++j) {
		for (int j = 0; j < NumberOfLStrip; ++j) { // back strip are the longitudinal one
			//StripCenter = Strip_1_1 + StripPitch * (i * U + j * V);
			StripCenter  = Strip_1_1 + StripLPitch*( j*U ) + StripTPitch*( i*V ); // or: Strip_1_1 + StripTPitch*( j*U ) + StripLPitch*( i*V  ) !!! To be checked as above!!!
			//StripCenter  = Strip_1_1 + StripTPitch*( j*U ) + StripLPitch*( i*V  ); // !!! To be checked as above!!!
			
			lineX.push_back(StripCenter.X());
			lineY.push_back(StripCenter.Y());
			lineZ.push_back(StripCenter.Z());
		}
		
		OneTelescopeStripPositionX.push_back(lineX);
		OneTelescopeStripPositionY.push_back(lineY);
		OneTelescopeStripPositionZ.push_back(lineZ);
	}
	m_StripPositionX.push_back(OneTelescopeStripPositionX);
	m_StripPositionY.push_back(OneTelescopeStripPositionY);
	m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
}


void TIssPhysics::AddTelescope(double theta, double phi, double distance,
							   double beta_u, double beta_v, double beta_w) {
	/*
	 //  m_NumberOfTelescope++;
	 m_NumberOfDetector++;
	 
	 double Pi = 3.141592654;
	 
	 // convert from degree to radian:
	 theta = theta * Pi / 180.;
	 phi   = phi * Pi / 180.;
	 
	 // Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that Y
	 // strip are allong X axis)
	 TVector3 U;
	 // Vector V on Telescope Face (parallele to X Strip)
	 TVector3 V;
	 // Vector W normal to Telescope Face (pointing CsI)
	 TVector3 W;
	 // Vector position of Telescope Face center
	 TVector3 C;
	 
	 C = TVector3(distance * sin(theta) * cos(phi),
	 distance * sin(theta) * sin(phi), distance * cos(theta));
	 
	 TVector3 P
	 = TVector3(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
	 
	 W = C.Unit();
	 U = W.Cross(P);
	 V = W.Cross(U);
	 
	 U = U.Unit();
	 V = V.Unit();
	 
	 U.Rotate(beta_u * Pi / 180., U);
	 V.Rotate(beta_u * Pi / 180., U);
	 
	 U.Rotate(beta_v * Pi / 180., V);
	 V.Rotate(beta_v * Pi / 180., V);
	 
	 U.Rotate(beta_w * Pi / 180., W);
	 V.Rotate(beta_w * Pi / 180., W);
	 
	 double Face          = 98; // mm
	 double NumberOfStrip = 128;
	 double StripPitch    = Face / NumberOfStrip; // mm
	 
	 vector<double> lineX;
	 vector<double> lineY;
	 vector<double> lineZ;
	 
	 vector<vector<double>> OneTelescopeStripPositionX;
	 vector<vector<double>> OneTelescopeStripPositionY;
	 vector<vector<double>> OneTelescopeStripPositionZ;
	 
	 double X, Y, Z;
	 
	 // Moving C to the 1.1 corner:
	 C.SetX(C.X() - (Face / 2 - StripPitch / 2) * (V.X() + U.X()));
	 C.SetY(C.Y() - (Face / 2 - StripPitch / 2) * (V.Y() + U.Y()));
	 C.SetZ(C.Z() - (Face / 2 - StripPitch / 2) * (V.Z() + U.Z()));
	 
	 for (int i = 0; i < 128; ++i) {
	 
	 lineX.clear();
	 lineY.clear();
	 lineZ.clear();
	 
	 for (int j = 0; j < 128; ++j) {
	 X = C.X() + StripPitch * (U.X() * i + V.X() * j);
	 Y = C.Y() + StripPitch * (U.Y() * i + V.Y() * j);
	 Z = C.Z() + StripPitch * (U.Z() * i + V.Z() * j);
	 
	 lineX.push_back(X);
	 lineY.push_back(Y);
	 lineZ.push_back(Z);
	 }
	 
	 OneTelescopeStripPositionX.push_back(lineX);
	 OneTelescopeStripPositionY.push_back(lineY);
	 OneTelescopeStripPositionZ.push_back(lineZ);
	 }
	 m_StripPositionX.push_back(OneTelescopeStripPositionX);
	 m_StripPositionY.push_back(OneTelescopeStripPositionY);
	 m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
	 */
	
}
/*
 TVector3 TIssPhysics::GetPositionOfInteraction(const int i) const {
 TVector3 Position
 = TVector3(GetStripPositionX(TelescopeNumber[i], Si_X[i], Si_Y[i]),
 GetStripPositionY(TelescopeNumber[i], Si_X[i], Si_Y[i]),
 GetStripPositionZ(TelescopeNumber[i], Si_X[i], Si_Y[i]));
 
 return (Position);
 }
 
 TVector3 TIssPhysics::GetTelescopeNormal(const int i) const {
 TVector3 U = TVector3(GetStripPositionX(TelescopeNumber[i], 128, 1),
 GetStripPositionY(TelescopeNumber[i], 128, 1),
 GetStripPositionZ(TelescopeNumber[i], 128, 1))
 
 - TVector3(GetStripPositionX(TelescopeNumber[i], 1, 1),
 GetStripPositionY(TelescopeNumber[i], 1, 1),
 GetStripPositionZ(TelescopeNumber[i], 1, 1));
 
 TVector3 V = TVector3(GetStripPositionX(TelescopeNumber[i], 128, 128),
 GetStripPositionY(TelescopeNumber[i], 128, 128),
 GetStripPositionZ(TelescopeNumber[i], 128, 128))
 
 - TVector3(GetStripPositionX(TelescopeNumber[i], 128, 1),
 GetStripPositionY(TelescopeNumber[i], 128, 1),
 GetStripPositionZ(TelescopeNumber[i], 128, 1));
 
 TVector3 Normal = U.Cross(V);
 
 return (Normal.Unit());
 }
 
 */






// Box
void TIssPhysics::AddBoxDetector(double Z){
	// BOX //
	double BOX_PCB_Width  = 26;  //61.10;
	double BOX_PCB_Length = 129; //104.00;
	double BOX_PCB_Thickness = 1;// 3.4;
	double BOX_PCB_Border_LongSide = 1;
	double BOX_PCB_Border_ShortSide = 2;
	
	// Single stage box case (DSSD only)
	double BOX_PCB_Slot_Width1 = BOX_PCB_Thickness;
	double BOX_PCB_Slot_Border1 = 4;
	double BOX_PCB_Slot_Deepness1 = BOX_PCB_Border_ShortSide;
	
	// BOX Wafer
	double BOX_ActiveWafer_Width  = 22; //48;
	double BOX_ActiveWafer_Length = 125;//72;
	double BOX_Wafer_Width  = 23;       // 52.20;
	double BOX_Wafer_Length = 126;      //76.20;
	
	int    BOX_Wafer_Front_NumberOfStrip = 128; //24 ;
	int    BOX_Wafer_Back_NumberOfStrip = 11; //22; //48 ;
	
	// Compute
	double BOX_LeftOver1 =  BOX_PCB_Length - BOX_PCB_Border_ShortSide - BOX_Wafer_Length - BOX_PCB_Slot_Border1 - BOX_PCB_Slot_Width1 ;
	double BOX_Exposed_Length1 = BOX_Wafer_Length + BOX_PCB_Slot_Border1 ;
	
	double BOX_CenterOffset1 = - 0.5 * BOX_PCB_Length+BOX_PCB_Border_ShortSide+0.5*BOX_Exposed_Length1;
	double BOX_DetectorSpacing1 = 0.5*BOX_Exposed_Length1+0.5*BOX_PCB_Slot_Width1;
	
	double BOX_Wafer_Width_Offset1 = -0.5*BOX_PCB_Width + BOX_PCB_Border_LongSide + 0.5*BOX_Wafer_Width;
	double BOX_Wafer_Length_Offset1 = -0.5*BOX_PCB_Length + BOX_PCB_Border_ShortSide + 0.5*BOX_Wafer_Length;
	
	double BOX_PCB_Slot_Position1 = 0.5*BOX_PCB_Length-BOX_LeftOver1 - 0.5*BOX_PCB_Slot_Width1;
	
	double StripPitchFront = BOX_ActiveWafer_Length/BOX_Wafer_Front_NumberOfStrip ; //mm
	double StripPitchBack  = BOX_ActiveWafer_Width/BOX_Wafer_Back_NumberOfStrip ; //mm
	m_BoxPitchBack = StripPitchBack;
	m_BoxPitchFront = StripPitchFront;
	
	// Double stage box case (DSSD+PAD) (the wafer is the same but the slot is different to accomodate the additional PAD)
	double PAD_PCB_Thickness = 3.4;
	
	double BOX_PCB_Slot_Width2 = BOX_PCB_Thickness + PAD_PCB_Thickness ;
	double BOX_PCB_Slot_Border2 = 2.7;
	double BOX_PCB_Slot_Deepness2 = BOX_PCB_Border_ShortSide;
	
	double BOX_LeftOver2 =  BOX_PCB_Length - BOX_PCB_Border_ShortSide - BOX_Wafer_Length - BOX_PCB_Slot_Border2 - BOX_PCB_Slot_Width2;
	double BOX_Exposed_Length2 = BOX_Wafer_Length + BOX_PCB_Slot_Border2 ;
	
	double BOX_CenterOffset2 = - 0.5*BOX_PCB_Length+BOX_PCB_Border_ShortSide + 0.5*BOX_Exposed_Length2;
	double BOX_DetectorSpacing2 = 0.5*BOX_Exposed_Length2 + 0.5*BOX_PCB_Thickness;
	
	double BOX_Wafer_Width_Offset2 = - 0.5*BOX_PCB_Width + BOX_PCB_Border_LongSide + 0.5*BOX_Wafer_Width;
	double BOX_Wafer_Length_Offset2 = - 0.5*BOX_PCB_Length + BOX_PCB_Border_ShortSide + 0.5*BOX_Wafer_Length;
	
	double BOX_PCB_Slot_Position2 = 0.5*BOX_PCB_Length-BOX_LeftOver2 - 0.5*BOX_PCB_Slot_Width2;
	
	
	double A1 = BOX_Exposed_Length1*0.5 -BOX_PCB_Slot_Border1- 0.5*StripPitchFront ;
	double B1 = BOX_DetectorSpacing1 - 0.5*BOX_PCB_Thickness;
	double Z1 = Z - BOX_Wafer_Width*0.5 + StripPitchBack*0.5 ;
	
	double A2 = BOX_Exposed_Length2*0.5 -BOX_PCB_Slot_Border2- 0.5*StripPitchFront ;
	double B2 = BOX_DetectorSpacing2 - 0.5*BOX_PCB_Thickness;
	double Z2 = Z + BOX_Wafer_Width*0.5 - StripPitchBack*0.5 ;
	
	TVector3 U; TVector3 V;TVector3 Strip_1_1;
	
	
	
	//cout << "cava" << endl;
	
	
	// To do : adapted for 24 detecteor rather than 4 and for 6 sides/orientation along z rather than 4
	for(int i = 0 ; i < 4 ; i++){
		m_NumberOfDetector++;
		if(Z<0){// Up Stream
			if(i==0)      {U=TVector3(1,0,0);V=TVector3(0,0,1);  Strip_1_1=TVector3( -A1 , B1  ,Z1); m_DetectorNormal.push_back(TVector3(0,-1,0));}
			else if(i==1) {U=TVector3(0,1,0);V=TVector3(0,0,1);  Strip_1_1=TVector3( -B1 , -A1 ,Z1); m_DetectorNormal.push_back(TVector3(1,0,0)) ;}
			else if(i==2) {U=TVector3(-1,0,0);V=TVector3(0,0,1); Strip_1_1=TVector3( A1  , -B1 ,Z1); m_DetectorNormal.push_back(TVector3(0,1,0)) ;}
			else if(i==3) {U=TVector3(0,-1,0);V=TVector3(0,0,1); Strip_1_1=TVector3( B1  , A1  ,Z1); m_DetectorNormal.push_back(TVector3(-1,0,0));}
		}
		
		else if(Z>0){//Down Stream
			if(i==0)      {U=TVector3(-1,0,0);V=TVector3(0,0,-1); Strip_1_1=TVector3( A2  ,B2  ,Z2); m_DetectorNormal.push_back(TVector3(0,-1,0));}
			else if(i==1) {U=TVector3(0,-1,0);V=TVector3(0,0,-1); Strip_1_1=TVector3( -B2 ,A2  ,Z2); m_DetectorNormal.push_back(TVector3(1,0,0)) ;}
			else if(i==2) {U=TVector3(1,0,0);V=TVector3(0,0,-1);  Strip_1_1=TVector3( -A2 ,-B2 ,Z2); m_DetectorNormal.push_back(TVector3(0,1,0)) ;}
			else if(i==3) {U=TVector3(0,1,0);V=TVector3(0,0,-1);  Strip_1_1=TVector3( B2  ,-A2 ,Z2); m_DetectorNormal.push_back(TVector3(-1,0,0));}
		}
		
		m_U.push_back(U);
		m_V.push_back(V);
		
		//cout << "cava" << endl;
		
		
		//   Buffer object to fill Position Array
		vector<double> lineX ; vector<double> lineY ; vector<double> lineZ ;
		
		vector< vector< double > >   OneBoxStripPositionX   ;
		vector< vector< double > >   OneBoxStripPositionY   ;
		vector< vector< double > >   OneBoxStripPositionZ   ;
		
		TVector3 StripCenter = Strip_1_1;
		for(int f = 0 ; f < BOX_Wafer_Front_NumberOfStrip ; f++){
			lineX.clear()   ;
			lineY.clear()   ;
			lineZ.clear()   ;
			
			for(int b = 0 ; b < BOX_Wafer_Back_NumberOfStrip ; b++){
				StripCenter = Strip_1_1 + ( StripPitchFront*f*U + StripPitchBack*b*V  );
				
				lineX.push_back( StripCenter.X() );
				lineY.push_back( StripCenter.Y() );
				lineZ.push_back( StripCenter.Z() );
			}
			
			OneBoxStripPositionX.push_back(lineX);
			OneBoxStripPositionY.push_back(lineY);
			OneBoxStripPositionZ.push_back(lineZ);
		}
		m_StripPositionX.push_back( OneBoxStripPositionX ) ;
		m_StripPositionY.push_back( OneBoxStripPositionY ) ;
		m_StripPositionZ.push_back( OneBoxStripPositionZ ) ;
	}
}


////////////////////////////////////////////////////////////////////////////////
TVector3 TIssPhysics::GetDetectorNormal( const int& i) const{
	return (m_DetectorNormal[DetectorNumber[i]-1]);
	
}
////////////////////////////////////////////////////////////////////////////////
TVector3 TIssPhysics::GetPositionOfInteraction(const int& i,bool random) const{
	static TVector3 Position ;
	
	//cout << "cava Pos" << endl;
	
	
	Position = TVector3 (  GetStripPositionX( DetectorNumber[i], Strip_Front[i], Strip_Back[i] ),
						 GetStripPositionY( DetectorNumber[i] , Strip_Front[i], Strip_Back[i] ),
						 GetStripPositionZ( DetectorNumber[i] , Strip_Front[i], Strip_Back[i] )) ;
	
	
	if(random){
		// Box Detector
		
		//cout << "DetNum: " << DetectorNumber[i] << endl;
		if(m_U[ DetectorNumber[i]-1].Mag()!=0){
			Position += m_V[ DetectorNumber[i]-1]*m_Rand->Uniform(-1,1)*m_BoxPitchBack*0.5;
			Position += m_U[ DetectorNumber[i]-1]*m_Rand->Uniform(-1,1)*m_BoxPitchFront*0.5;
		}
		
	}
	
	return Position ;
	
}
////////////////////////////////////////////////////////////////////////////////
double TIssPhysics::GetDeadLayer(const int& i ) const{
	return DeadLayer[i];
}
////////////////////////////////////////////////////////////////////////////////
void TIssPhysics::InitializeStandardParameter()
{
	//   Enable all channel
	vector< bool > ChannelStatus;
	m_FrontChannelStatus.clear()    ;
	m_BackChannelStatus.clear()    ;
	
	
	ChannelStatus.resize(128,true);
	for(int i = 0 ; i < m_NumberOfDetector ; i++){
		m_FrontChannelStatus[i] = ChannelStatus;
	}
	
	ChannelStatus.resize(11,true);
	for(int i = 0 ; i < m_NumberOfDetector ; i++){
		m_BackChannelStatus[i] = ChannelStatus;
	}
	
	
	m_MaximumStripMultiplicityAllowed = m_NumberOfDetector   ;
	
	return;
}


///////////////////////////////////////////////////////////////////////////
namespace ISS_LOCAL{
//   DSSD
//   Front
double fStrip_Front_E(const TIssData* m_EventData , const int& i){
	static CalibrationManager* Cal = CalibrationManager::getInstance();
	static string name ;
	name = "ISS/D" + NPL::itoa( m_EventData->GetFront_DetectorNbr(i) ) + "_STRIP_FRONT" + NPL::itoa( m_EventData->GetFront_StripNbr(i) ) + "_E";
	return Cal->ApplyCalibration(name,m_EventData->GetFront_Energy(i) );
}

double fStrip_Front_T(const TIssData* m_EventData , const int& i){
	static CalibrationManager* Cal = CalibrationManager::getInstance();
	static string name ;
	name ="ISS/D" + NPL::itoa( m_EventData->GetFront_DetectorNbr(i) ) + "_STRIP_FRONT" + NPL::itoa( m_EventData->GetFront_StripNbr(i) ) +"_T";
	
	return Cal->ApplyCalibration(name, m_EventData->GetFront_TimeCFD(i) );
}

//   Back
double fStrip_Back_E(const TIssData* m_EventData , const int& i){
	static CalibrationManager* Cal = CalibrationManager::getInstance();
	static string name ;
	name =  "ISS/D" + NPL::itoa( m_EventData->GetBack_DetectorNbr(i) ) + "_STRIP_BACK" + NPL::itoa( m_EventData->GetBack_StripNbr(i)) +"_E";
	
	return Cal->ApplyCalibration(name, m_EventData->GetBack_Energy(i) );
}

double fStrip_Back_T(const TIssData* m_EventData , const int& i){
	static CalibrationManager* Cal = CalibrationManager::getInstance();
	static string name ;
	name = "ISS/D" + NPL::itoa( m_EventData->GetBack_DetectorNbr(i) ) + "_STRIP_BACK" + NPL::itoa( m_EventData->GetBack_StripNbr(i) ) +"_T";
	
	return Cal->ApplyCalibration(name, m_EventData->GetFront_TimeCFD(i));
}


}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TIssPhysics::Construct(){
	return (NPL::VDetector*) new TIssPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_iss{
public:
	proxy_iss(){
		NPL::DetectorFactory::getInstance()->AddToken("ISSDet","Iss");
		NPL::DetectorFactory::getInstance()->AddDetector("ISSDet",TIssPhysics::Construct);
	}
};

proxy_iss p_iss;
}




/* from MUST2:
 ////////////////////////////////////////////////////////////////////////////
 int TIssPhysics::CheckEvent(int N) {
 // Bad event
 if (m_NMatchDet[N] > m_StripXMultDet[N]
 || m_NMatchDet[N] > m_StripYMultDet[N]) {
 return  2;
 }
 // Good event
 else {
 return 1;
 }
 }
 
 ////////////////////////////////////////////////////////////////////////////
 bool TIssPhysics::IsValidChannel(const int& DetectorType,
 const int& telescope, const int& channel) {
 if (DetectorType == 0)
 return *(m_XChannelStatus[telescope - 1].begin() + channel - 1);
 
 else if (DetectorType == 1)
 return *(m_YChannelStatus[telescope - 1].begin() + channel - 1);
 
 else if (DetectorType == 2)
 return *(m_SiLiChannelStatus[telescope - 1].begin() + channel - 1);
 
 else if (DetectorType == 3)
 return *(m_CsIChannelStatus[telescope - 1].begin() + channel - 1);
 
 else
 return false;
 }
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::ReadAnalysisConfig() {
 bool ReadingStatus = false;
 
 // path to file
 string FileName = "./configs/ConfigIss.dat";
 
 // open analysis config file
 ifstream AnalysisConfigFile;
 AnalysisConfigFile.open(FileName.c_str());
 
 if (!AnalysisConfigFile.is_open()) {
 cout << " No ConfigIss.dat found: Default parameters loaded for "
 "Analysis "
 << FileName << endl;
 return;
 }
 cout << " Loading user parameters for Analysis from ConfigIss.dat " << endl;
 
 // Save it in a TAsciiFile
 TAsciiFile* asciiConfig
 = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
 asciiConfig->AppendLine("%%% ConfigIss.dat %%%");
 asciiConfig->Append(FileName.c_str());
 asciiConfig->AppendLine("");
 // read analysis config file
 string LineBuffer, DataBuffer, whatToDo;
 while (!AnalysisConfigFile.eof()) {
 // Pick-up next line
 getline(AnalysisConfigFile, LineBuffer);
 
 // search for "header"
 if (LineBuffer.compare(0, 11, "ConfigIss") == 0)
 ReadingStatus = true;
 
 // loop on tokens and data
 while (ReadingStatus) {
 
 whatToDo = "";
 AnalysisConfigFile >> whatToDo;
 
 // Search for comment symbol (%)
 if (whatToDo.compare(0, 1, "%") == 0) {
 AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
 }
 
 else if (whatToDo == "MAX_STRIP_MULTIPLICITY") {
 AnalysisConfigFile >> DataBuffer;
 m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str());
 cout << "MAXIMUN STRIP MULTIPLICITY "
 << m_MaximumStripMultiplicityAllowed << endl;
 }
 
 else if (whatToDo == "STRIP_ENERGY_MATCHING_SIGMA") {
 AnalysisConfigFile >> DataBuffer;
 m_StripEnergyMatchingSigma = atof(DataBuffer.c_str());
 cout << "STRIP ENERGY MATCHING SIGMA " << m_StripEnergyMatchingSigma
 << endl;
 }
 
 else if (whatToDo == "STRIP_ENERGY_MATCHING_NUMBER_OF_SIGMA") {
 AnalysisConfigFile >> DataBuffer;
 m_StripEnergyMatchingNumberOfSigma = atoi(DataBuffer.c_str());
 cout << "STRIP ENERGY MATCHING NUMBER OF SIGMA "
 << m_StripEnergyMatchingNumberOfSigma << endl;
 }
 
 else if (whatToDo == "DISABLE_ALL") {
 AnalysisConfigFile >> DataBuffer;
 cout << whatToDo << "  " << DataBuffer << endl;
 int          telescope = atoi(DataBuffer.substr(2, 1).c_str());
 vector<bool> ChannelStatus;
 ChannelStatus.resize(128, false);
 m_XChannelStatus[telescope - 1] = ChannelStatus;
 m_YChannelStatus[telescope - 1] = ChannelStatus;
 ChannelStatus.resize(16, false);
 m_SiLiChannelStatus[telescope - 1] = ChannelStatus;
 m_CsIChannelStatus[telescope - 1]  = ChannelStatus;
 }
 
 else if (whatToDo == "DISABLE_CHANNEL") {
 AnalysisConfigFile >> DataBuffer;
 cout << whatToDo << "  " << DataBuffer << endl;
 int telescope = atoi(DataBuffer.substr(2, 1).c_str());
 int channel   = -1;
 if (DataBuffer.compare(3, 4, "STRX") == 0) {
 channel = atoi(DataBuffer.substr(7).c_str());
 *(m_XChannelStatus[telescope - 1].begin() + channel - 1) = false;
 }
 
 else if (DataBuffer.compare(3, 4, "STRY") == 0) {
 channel = atoi(DataBuffer.substr(7).c_str());
 *(m_YChannelStatus[telescope - 1].begin() + channel - 1) = false;
 }
 
 
 else
 cout << "Warning: detector type for Iss unknown!" << endl;
 
 }
 
 else if (whatToDo == "TAKE_E_Y") {
 m_Take_E_Y = true;
 cout << whatToDo << endl;
 }
 
 else if (whatToDo == "TAKE_T_Y") {
 m_Take_T_Y = true;
 cout << whatToDo << endl;
 }
 
 else if (whatToDo == "TAKE_E_X") {
 m_Take_E_Y = false;
 cout << whatToDo << endl;
 }
 
 else if (whatToDo == "TAKE_T_X") {
 m_Take_T_Y = false;
 cout << whatToDo << endl;
 }
 
 else if (whatToDo == "IGNORE_NOT_MATCHING_CSI") {
 m_Ignore_not_matching_CsI = true;
 cout << whatToDo << endl;
 }
 
 
 else if (whatToDo == "SI_X_E_RAW_THRESHOLD") {
 AnalysisConfigFile >> DataBuffer;
 m_Si_X_E_RAW_Threshold = atof(DataBuffer.c_str());
 cout << whatToDo << " " << m_Si_X_E_RAW_Threshold << endl;
 }
 
 else if (whatToDo == "SI_Y_E_RAW_THRESHOLD") {
 AnalysisConfigFile >> DataBuffer;
 m_Si_Y_E_RAW_Threshold = atof(DataBuffer.c_str());
 cout << whatToDo << " " << m_Si_Y_E_RAW_Threshold << endl;
 }
 
 
 else if (whatToDo == "SI_X_E_THRESHOLD") {
 AnalysisConfigFile >> DataBuffer;
 m_Si_X_E_Threshold = atof(DataBuffer.c_str());
 cout << whatToDo << " " << m_Si_X_E_Threshold << endl;
 }
 
 else if (whatToDo == "SI_Y_E_THRESHOLD") {
 AnalysisConfigFile >> DataBuffer;
 m_Si_Y_E_Threshold = atof(DataBuffer.c_str());
 cout << whatToDo << " " << m_Si_Y_E_Threshold << endl;
 }
 
 
 else {
 ReadingStatus = false;
 }
 }
 }
 }
 
 ///////////////////////////////////////////////////////////////////////////
 bool TIssPhysics::Match_Si_SiLi(int X, int Y, int PadNbr) {
 
 // remove the central part and surrounding
 if (X < 8 || X > 120 || (Y < 68 && Y > 60))
 return false;
 
 if (abs(m_SiLi_MatchingX[PadNbr - 1] - X) < m_SiLi_Size / 2.
 && abs(m_SiLi_MatchingY[PadNbr - 1] - Y) < m_SiLi_Size / 2.)
 return true;
 
 else
 return false;
 }
 
 ///////////////////////////////////////////////////////////////////////////
 bool TIssPhysics::Match_Si_CsI(int X, int Y, int CristalNbr) {
 
 if (abs(m_CsI_MatchingX[CristalNbr - 1] - X) < (double)m_CsI_Size / 2.
 && abs(m_CsI_MatchingY[CristalNbr - 1] - Y) < (double)m_CsI_Size / 2.) {
 return true;
 }
 
 else
 return false;
 }
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::Clear() {
 EventMultiplicity = 0;
 
 m_match_type.clear();
 
 m_StripXMultDet.clear();
 m_StripYMultDet.clear();
 m_NMatchDet.clear();
 
 TelescopeNumber.clear();
 EventType.clear();
 TotalEnergy.clear();
 
 // Si X
 Si_E.clear();
 Si_T.clear();
 Si_X.clear();
 Si_Y.clear();
 
 
 Si_EX.clear();
 Si_TX.clear();
 Si_EY.clear();
 Si_TY.clear();
 TelescopeNumber_X.clear();
 TelescopeNumber_Y.clear();
 }
 ///////////////////////////////////////////////////////////////////////////
 
 void TIssPhysics::ReadCalibrationRun() {
 m_StripXEMult = m_EventData->GetMMStripXEMult();
 m_StripYEMult = m_EventData->GetMMStripYEMult();
 m_StripXTMult = m_EventData->GetMMStripXTMult();
 m_StripYTMult = m_EventData->GetMMStripYTMult();
 
 
 //   X
 //   E
 for (unsigned int i = 0; i < m_StripXEMult; ++i) {
 TelescopeNumber_X.push_back(m_EventData->GetMMStripXEDetectorNbr(i));
 Si_EX.push_back(fSi_X_E(m_EventData, i));
 Si_X.push_back(m_EventData->GetMMStripXEStripNbr(i));
 }
 //   T
 for (unsigned int i = 0; i < m_StripXTMult; ++i) {
 TelescopeNumber_X.push_back(m_EventData->GetMMStripXTDetectorNbr(i));
 Si_TX.push_back(fSi_X_T(m_EventData, i));
 Si_X.push_back(m_EventData->GetMMStripXTStripNbr(i));
 }
 
 //   Y
 //   E
 for (unsigned int i = 0; i < m_StripYEMult; ++i) {
 TelescopeNumber_Y.push_back(m_EventData->GetMMStripYEDetectorNbr(i));
 Si_EY.push_back(fSi_Y_E(m_EventData, i));
 Si_Y.push_back(m_EventData->GetMMStripYEStripNbr(i));
 }
 
 //   T
 for (unsigned int i = 0; i < m_StripYTMult; ++i) {
 TelescopeNumber_Y.push_back(m_EventData->GetMMStripYTDetectorNbr(i));
 Si_TY.push_back(fSi_Y_T(m_EventData, i));
 Si_Y.push_back(m_EventData->GetMMStripYTStripNbr(i));
 }
 }
 
 
 ////   Innherited from VDetector Class   ////
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::ReadConfiguration(NPL::InputParser parser) {
 vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ISSDet");
 if (NPOptionManager::getInstance()->GetVerboseLevel())
 cout << "//// " << blocks.size() << " Detector found " << endl;
 
 // Cartesian Case
 vector<string> cart
 = {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128", "SI"};
 // Spherical Case
 vector<string> sphe = {"R", "THETA", "PHI", "BETA", "SI"};
 
 for (unsigned int i = 0; i < blocks.size(); i++) {
 if (blocks[i]->HasTokenList(cart)) {
 if (NPOptionManager::getInstance()->GetVerboseLevel())
 cout << endl << "////  Iss Telescope " << i + 1 << endl;
 TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
 TVector3 B = blocks[i]->GetTVector3("X128_Y1", "mm");
 TVector3 C = blocks[i]->GetTVector3("X1_Y128", "mm");
 TVector3 D = blocks[i]->GetTVector3("X128_Y128", "mm");
 AddTelescope(A, B, C, D);
 }
 
 else if (blocks[i]->HasTokenList(sphe)) {
 if (NPOptionManager::getInstance()->GetVerboseLevel())
 cout << endl << "////  Iss Telescope " << i + 1 << endl;
 
 double         Theta = blocks[i]->GetDouble("THETA", "deg");
 double         Phi   = blocks[i]->GetDouble("PHI", "deg");
 double         R     = blocks[i]->GetDouble("R", "mm");
 vector<double> beta  = blocks[i]->GetVectorDouble("BETA", "deg");
 AddTelescope(Theta, Phi, R, beta[0], beta[1], beta[2]);
 
 }
 
 else {
 cout << "ERROR: Missing token for M2Telescope blocks, check your "
 "input "
 "file"
 << endl;
 exit(1);
 }
 }
 
 // m_MultEvt = new TIssMultTelescope[m_NumberOfTelescope];
 
 InitializeStandardParameter();
 ReadAnalysisConfig();
 }
 //////////////////////////////////////////////////////////////////////////
 void TIssPhysics::InitSpectra() {
 m_Spectra = new TIssSpectra(m_NumberOfTelescope);
 }
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::FillSpectra() {
 m_Spectra->FillRawSpectra(m_EventData);
 m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
 m_Spectra->FillPhysicsSpectra(m_EventPhysics);
 }
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::CheckSpectra() { m_Spectra->CheckSpectra(); }
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::ClearSpectra() {
 // To be done
 }
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::WriteSpectra() {
 if (m_Spectra)
 m_Spectra->WriteSpectra();
 }
 
 ///////////////////////////////////////////////////////////////////////////
 map<string, TH1*> TIssPhysics::GetSpectra() {
 if (m_Spectra)
 return m_Spectra->GetMapHisto();
 else {
 map<string, TH1*> empty;
 return empty;
 }
 }
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::AddParameterToCalibrationManager() {
 CalibrationManager* Cal = CalibrationManager::getInstance();
 // Good for simulation, close to typical values
 vector<double> standardX    = {-63, 63. / 8192.};
 vector<double> standardY    = {63, -63. / 8192.};
 vector<double> standardT    = {-1000, 1000. / 8192.};
 
 for (int i = 0; i < m_NumberOfTelescope; ++i) {
 
 for (int j = 0; j < 128; ++j) {
 Cal->AddParameter(
 "Iss", "T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_E",
 "Iss_T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_E",
 standardX);
 Cal->AddParameter(
 "Iss", "T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_E",
 "Iss_T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_E",
 standardY);
 Cal->AddParameter(
 "Iss", "T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_T",
 "Iss_T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_T",
 standardT);
 Cal->AddParameter(
 "Iss", "T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_T",
 "Iss_T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_T",
 standardT);
 }
 
 }
 
 return;
 }
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::InitializeRootInputRaw() {
 TChain* inputChain = RootInput::getInstance()->GetChain();
 inputChain->SetBranchStatus("ISS", true);
 inputChain->SetBranchStatus("fMM_*", true);
 inputChain->SetBranchAddress("ISS", &m_EventData);
 }
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::InitializeRootInputPhysics() {
 TChain* inputChain = RootInput::getInstance()->GetChain();
 inputChain->SetBranchStatus("ISS", true);
 inputChain->SetBranchStatus("EventMultiplicity", true);
 inputChain->SetBranchStatus("EventType", true);
 inputChain->SetBranchStatus("DetNumber", true);
 inputChain->SetBranchStatus("Si_E", true);
 inputChain->SetBranchStatus("Si_T", true);
 inputChain->SetBranchStatus("Si_X", true);
 inputChain->SetBranchStatus("Si_Y", true);
 inputChain->SetBranchStatus("Si_EX", true);
 inputChain->SetBranchStatus("Si_TX", true);
 inputChain->SetBranchStatus("Si_EY", true);
 inputChain->SetBranchStatus("Si_TY", true);
 inputChain->SetBranchStatus("DetNumber_X", true);
 inputChain->SetBranchStatus("DetNumber_Y", true);
 inputChain->SetBranchStatus("TotalEnergy", true);
 inputChain->SetBranchAddress("ISS", &m_EventPhysics);
 }
 
 ///////////////////////////////////////////////////////////////////////////
 void TIssPhysics::InitializeRootOutput() {
 TTree* outputTree = RootOutput::getInstance()->GetTree();
 outputTree->Branch("ISS", "TIssPhysics", &m_EventPhysics);
 }
 
 /////   Specific to IssArray   ////
 
 void TIssPhysics::AddTelescope(TVector3 C_X1_Y1, TVector3 C_X128_Y1,
 TVector3 C_X1_Y128, TVector3 C_X128_Y128) {
 // To avoid warning
 C_X128_Y128 *= 1;
 
 m_NumberOfTelescope++;
 
 //   Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that
 //   Y strip are allong X axis)
 TVector3 U      = C_X128_Y1 - C_X1_Y1;
 double   Ushift = (U.Mag() - 98) / 2.;
 U               = U.Unit();
 //   Vector V on Telescope Face (parallele to X Strip)
 TVector3 V      = C_X1_Y128 - C_X1_Y1;
 double   Vshift = (V.Mag() - 98) / 2.;
 V               = V.Unit();
 
 //   Position Vector of Strip Center
 TVector3 StripCenter = TVector3(0, 0, 0);
 //   Position Vector of X=1 Y=1 Strip
 TVector3 Strip_1_1;
 
 //   Geometry Parameter
 double Face          = 98; // mm
 double NumberOfStrip = 128;
 double StripPitch    = Face / NumberOfStrip; // mm
 //   Buffer object to fill Position Array
 vector<double> lineX;
 vector<double> lineY;
 vector<double> lineZ;
 
 vector<vector<double>> OneTelescopeStripPositionX;
 vector<vector<double>> OneTelescopeStripPositionY;
 vector<vector<double>> OneTelescopeStripPositionZ;
 
 //   Moving StripCenter to 1.1 corner:
 Strip_1_1 = C_X1_Y1 + (U + V) * (StripPitch / 2.);
 Strip_1_1 += U * Ushift + V * Vshift;
 
 for (int i = 0; i < 128; ++i) {
 lineX.clear();
 lineY.clear();
 lineZ.clear();
 
 for (int j = 0; j < 128; ++j) {
 StripCenter = Strip_1_1 + StripPitch * (i * U + j * V);
 lineX.push_back(StripCenter.X());
 lineY.push_back(StripCenter.Y());
 lineZ.push_back(StripCenter.Z());
 }
 
 OneTelescopeStripPositionX.push_back(lineX);
 OneTelescopeStripPositionY.push_back(lineY);
 OneTelescopeStripPositionZ.push_back(lineZ);
 }
 m_StripPositionX.push_back(OneTelescopeStripPositionX);
 m_StripPositionY.push_back(OneTelescopeStripPositionY);
 m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
 }
 
 void TIssPhysics::InitializeStandardParameter() {
 //   Enable all channel
 vector<bool> ChannelStatus;
 m_XChannelStatus.clear();
 m_YChannelStatus.clear();
 
 
 ChannelStatus.resize(128, true);
 for (int i = 0; i < m_NumberOfTelescope; ++i) {
 m_XChannelStatus[i] = ChannelStatus;
 m_YChannelStatus[i] = ChannelStatus;
 }
 
 
 m_MaximumStripMultiplicityAllowed = m_NumberOfTelescope;
 
 return;
 }
 
 void TIssPhysics::AddTelescope(double theta, double phi, double distance,
 double beta_u, double beta_v, double beta_w) {
 
 m_NumberOfTelescope++;
 
 double Pi = 3.141592654;
 
 // convert from degree to radian:
 theta = theta * Pi / 180.;
 phi   = phi * Pi / 180.;
 
 // Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that Y
 // strip are allong X axis)
 TVector3 U;
 // Vector V on Telescope Face (parallele to X Strip)
 TVector3 V;
 // Vector W normal to Telescope Face (pointing CsI)
 TVector3 W;
 // Vector position of Telescope Face center
 TVector3 C;
 
 C = TVector3(distance * sin(theta) * cos(phi),
 distance * sin(theta) * sin(phi), distance * cos(theta));
 
 TVector3 P
 = TVector3(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
 
 W = C.Unit();
 U = W.Cross(P);
 V = W.Cross(U);
 
 U = U.Unit();
 V = V.Unit();
 
 U.Rotate(beta_u * Pi / 180., U);
 V.Rotate(beta_u * Pi / 180., U);
 
 U.Rotate(beta_v * Pi / 180., V);
 V.Rotate(beta_v * Pi / 180., V);
 
 U.Rotate(beta_w * Pi / 180., W);
 V.Rotate(beta_w * Pi / 180., W);
 
 double Face          = 98; // mm
 double NumberOfStrip = 128;
 double StripPitch    = Face / NumberOfStrip; // mm
 
 vector<double> lineX;
 vector<double> lineY;
 vector<double> lineZ;
 
 vector<vector<double>> OneTelescopeStripPositionX;
 vector<vector<double>> OneTelescopeStripPositionY;
 vector<vector<double>> OneTelescopeStripPositionZ;
 
 double X, Y, Z;
 
 // Moving C to the 1.1 corner:
 C.SetX(C.X() - (Face / 2 - StripPitch / 2) * (V.X() + U.X()));
 C.SetY(C.Y() - (Face / 2 - StripPitch / 2) * (V.Y() + U.Y()));
 C.SetZ(C.Z() - (Face / 2 - StripPitch / 2) * (V.Z() + U.Z()));
 
 for (int i = 0; i < 128; ++i) {
 
 lineX.clear();
 lineY.clear();
 lineZ.clear();
 
 for (int j = 0; j < 128; ++j) {
 X = C.X() + StripPitch * (U.X() * i + V.X() * j);
 Y = C.Y() + StripPitch * (U.Y() * i + V.Y() * j);
 Z = C.Z() + StripPitch * (U.Z() * i + V.Z() * j);
 
 lineX.push_back(X);
 lineY.push_back(Y);
 lineZ.push_back(Z);
 }
 
 OneTelescopeStripPositionX.push_back(lineX);
 OneTelescopeStripPositionY.push_back(lineY);
 OneTelescopeStripPositionZ.push_back(lineZ);
 }
 m_StripPositionX.push_back(OneTelescopeStripPositionX);
 m_StripPositionY.push_back(OneTelescopeStripPositionY);
 m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
 }
 
 TVector3 TIssPhysics::GetPositionOfInteraction(const int i) const {
 TVector3 Position
 = TVector3(GetStripPositionX(TelescopeNumber[i], Si_X[i], Si_Y[i]),
 GetStripPositionY(TelescopeNumber[i], Si_X[i], Si_Y[i]),
 GetStripPositionZ(TelescopeNumber[i], Si_X[i], Si_Y[i]));
 
 return (Position);
 }
 
 TVector3 TIssPhysics::GetTelescopeNormal(const int i) const {
 TVector3 U = TVector3(GetStripPositionX(TelescopeNumber[i], 128, 1),
 GetStripPositionY(TelescopeNumber[i], 128, 1),
 GetStripPositionZ(TelescopeNumber[i], 128, 1))
 
 - TVector3(GetStripPositionX(TelescopeNumber[i], 1, 1),
 GetStripPositionY(TelescopeNumber[i], 1, 1),
 GetStripPositionZ(TelescopeNumber[i], 1, 1));
 
 TVector3 V = TVector3(GetStripPositionX(TelescopeNumber[i], 128, 128),
 GetStripPositionY(TelescopeNumber[i], 128, 128),
 GetStripPositionZ(TelescopeNumber[i], 128, 128))
 
 - TVector3(GetStripPositionX(TelescopeNumber[i], 128, 1),
 GetStripPositionY(TelescopeNumber[i], 128, 1),
 GetStripPositionZ(TelescopeNumber[i], 128, 1));
 
 TVector3 Normal = U.Cross(V);
 
 return (Normal.Unit());
 }
 
 ///////////////////////////////////////////////////////////////////////////
 namespace ISS_LOCAL {
 //   DSSD
 //   X
 double fSi_X_E(const TIssData* m_EventData, const int& i) {
 static string name;
 name = "ISS/T";
 name += NPL::itoa(m_EventData->GetMMStripXEDetectorNbr(i));
 name += "_Si_X";
 name += NPL::itoa(m_EventData->GetMMStripXEStripNbr(i));
 name += "_E";
 return CalibrationManager::getInstance()->ApplyCalibration(
 name, m_EventData->GetMMStripXEEnergy(i));
 }
 
 double fSi_X_T(const TIssData* m_EventData, const int& i) {
 static string name;
 name = "ISS/T";
 name += NPL::itoa(m_EventData->GetMMStripXTDetectorNbr(i));
 name += "_Si_X";
 name += NPL::itoa(m_EventData->GetMMStripXTStripNbr(i));
 name += "_T";
 return CalibrationManager::getInstance()->ApplyCalibration(
 name, m_EventData->GetMMStripXTTime(i));
 }
 
 //   Y
 double fSi_Y_E(const TIssData* m_EventData, const int& i) {
 static string name;
 name = "ISS/T";
 name += NPL::itoa(m_EventData->GetMMStripYEDetectorNbr(i));
 name += "_Si_Y";
 name += NPL::itoa(m_EventData->GetMMStripYEStripNbr(i));
 name += "_E";
 return CalibrationManager::getInstance()->ApplyCalibration(
 name, m_EventData->GetMMStripYEEnergy(i));
 }
 
 double fSi_Y_T(const TIssData* m_EventData, const int& i) {
 static string name;
 name = "ISS/T";
 name += NPL::itoa(m_EventData->GetMMStripYTDetectorNbr(i));
 name += "_Si_Y";
 name += NPL::itoa(m_EventData->GetMMStripYTStripNbr(i));
 name += "_T";
 return CalibrationManager::getInstance()->ApplyCalibration(
 name, m_EventData->GetMMStripYTTime(i));
 }
 
 } // namespace ISS_LOCAL
 
 ////////////////////////////////////////////////////////////////////////////////
 //            Construct Method to be pass to the DetectorFactory //
 ////////////////////////////////////////////////////////////////////////////////
 NPL::VDetector* TIssPhysics::Construct() {
 return (NPL::VDetector*)new TIssPhysics();
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 //            Registering the construct method to the factory //
 ////////////////////////////////////////////////////////////////////////////////
 extern "C" {
 class proxy_iss {
 public:
 proxy_iss() {
 NPL::DetectorFactory::getInstance()->AddToken("ISSDet", "ISS");
 NPL::DetectorFactory::getInstance()->AddDetector("ISSDet",
 TIssPhysics::Construct);
 }
 };
 
 proxy_iss p_iss;
 }
 */
