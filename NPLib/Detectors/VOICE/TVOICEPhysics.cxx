/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sunghan Bae  contact address: shbae2703@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold VOICE Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TVOICEPhysics.h"

//   STL
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"

//   ROOT
#include "TChain.h"
#include "NPInputParser.h"
#include "NPSystemOfUnits.h"
#include "TRotation.h"
#include "TAsciiFile.h"

using namespace NPUNITS;
using namespace std;

ClassImp(TVOICEPhysics)

///////////////////////////////////////////////////////////////////////////
TVOICEPhysics::TVOICEPhysics()
   : m_EventData(new TVOICEData),
     m_PreTreatedData(new TVOICEData),
     m_EventPhysics(this){
//     m_Spectra(0),
//   m_E_RAW_Threshold(0), // adc channels
// m_E_Threshold(0),     // MeV
//     m_NumberOfDetectors(0) {
}

TVOICEPhysics::~TVOICEPhysics() {}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TVOICEPhysics::AddDetector(string Type, double zpos, double Tilt, double Rot, double radi, double pitch){
	m_Type.push_back(Type);
	m_zpos.push_back(zpos);	
	m_Tilt.push_back(Tilt);
	m_Rot.push_back(Rot);	
	m_radi.push_back(radi);	
	m_pitch.push_back(pitch);		
} 

///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::BuildPhysicalEvent() {
  Clear();
	nhit = m_EventData->GetMult();
	for(Int_t i= 0 ; i < nhit ; i++){
		type[i] = m_EventData->GetType(i);
		E[i] = m_EventData->GetEnergy(i);
		X[i] = m_EventData->GetX(i);
		Y[i] = m_EventData->GetY(i);
		Z[i] = m_EventData->GetZ(i);
		Time[i] = m_EventData->GetTime(i);
		Theta[i] = m_EventData->GetTheta(i);
		Phi[i] = m_EventData->GetPhi(i);
		Atom_Z[i] = m_EventData->GetAtom_Z(i);
		A[i] = m_EventData->GetA(i);
		E_f[i] = m_EventData->GetEnergy_f(i);
		X_f[i] = m_EventData->GetX_f(i);
		Y_f[i] = m_EventData->GetY_f(i);
		Z_f[i] = m_EventData->GetZ_f(i);
		Time_f[i] = m_EventData->GetTime_f(i);
		Theta_f[i] = m_EventData->GetTheta_f(i);
		Phi_f[i] = m_EventData->GetPhi_f(i);
	}
}

void TVOICEPhysics::Clear(){
	nhit=0;
}

///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::ReadConfiguration(NPL::InputParser parser) {
	vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("VOICE");
	if (NPOptionManager::getInstance()->GetVerboseLevel())
		cout<<"//// "<< blocks.size() <<" detectors found "<<endl;

	bool FG_flag = false;
	bool Cathode_flag = false;
	vector<string> voice_vec = {"Type"};
	vector<string> voice_gas_vec = {"gas","pressure"};
	for(unsigned int i = 0 ; i < blocks.size() ; i++){

		int count_gas_info =0;
		if(blocks[i]->HasTokenList(voice_vec)){
			if(NPOptionManager::getInstance()->GetVerboseLevel())
				cout<< endl <<"\\\\ VOICE " << i+1 <<endl;

			string Type = blocks[i]->GetString("Type");
			double zpos = -9999;
			if(blocks[i]->HasToken("z-pos")) zpos = blocks[i]->GetDouble("z-pos","mm");
			double radi = 0;
			if(blocks[i]->HasToken("radi"))radi = blocks[i]->GetDouble("radi","micrometer");
			double pitch = 0;
			if(blocks[i]->HasToken("pitch"))pitch = blocks[i]->GetDouble("pitch","mm");
			double Tilt = 0;
			if(blocks[i]->HasToken("Tilt"))Tilt = blocks[i]->GetDouble("Tilt","deg");
			double zRot = 0;
			if(blocks[i]->HasToken("zRot"))zRot = blocks[i]->GetDouble("zRot","deg");

			if(blocks[i]->HasToken("AFG_gap")){
				if(!FG_flag) {
					m_AFG_gap = blocks[i]->GetDouble("AFG_gap","mm");
					FG_flag=true;
				}else{
					cout<<"============= CHECK!!!!!!!!! Multiple FG information input =============="<<endl;

				}	
			}

			if(blocks[i]->HasToken("AC_gap")){
				if(!Cathode_flag) {
					m_AC_gap = blocks[i]->GetDouble("AC_gap","mm");
					Cathode_flag=true;
				}else{
					cout<<"============= CHECK!!!!!!!!! Multiple FG information input =============="<<endl;
				}
			}

			AddDetector(Type,zpos,Tilt,zRot,radi,pitch);

		}else if (blocks[i]->HasTokenList(voice_gas_vec)){
			if(count_gas_info>0) {cout<<"ERROR: Multiple gas info."<<endl; exit(1);}
			m_gas = blocks[i]->GetString("gas");
			m_pressure = blocks[i]->GetDouble("pressure","atmosphere");
			count_gas_info++;
		}else{
			cout << "ERROR: check your input file formatting " << endl;
			exit(1);
		}
	}
	cout<< "read complete" <<endl;
}



///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::InitSpectra() {
	//	m_Spectra = new TVOICESpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::FillSpectra() {
	//	m_Spectra -> FillRawSpectra(m_EventData);
	//	m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
	//	m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::CheckSpectra() {
	//	m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::ClearSpectra() {
	// To be done
}



///////////////////////////////////////////////////////////////////////////
//map< string , TH1*> TVOICEPhysics::GetSpectra() {
//	if(m_Spectra)
//		return m_Spectra->GetMapHisto();
//	else{
//		map< string , TH1*> empty;
//		return empty;
//	}
//}

///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::WriteSpectra() {
	//	m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::AddParameterToCalibrationManager() {
	//	CalibrationManager* Cal = CalibrationManager::getInstance();
	//	for (int i = 0; i < m_NumberOfDetectors; ++i) {
	//		Cal->AddParameter("VOICE", "D"+ NPL::itoa(i+1)+"_ENERGY","VOICE_D"+ NPL::itoa(i+1)+"_ENERGY");
	//		Cal->AddParameter("VOICE", "D"+ NPL::itoa(i+1)+"_TIME","VOICE_D"+ NPL::itoa(i+1)+"_TIME");
	//	}
}



///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::InitializeRootInputRaw() {
	TChain* inputChain = RootInput::getInstance()->GetChain();
	inputChain->SetBranchStatus("VOICE",  true );
	inputChain->SetBranchAddress("VOICE", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::InitializeRootInputPhysics() {
	TChain* inputChain = RootInput::getInstance()->GetChain();
	inputChain->SetBranchStatus("VOICE",  true );
	inputChain->SetBranchAddress("VOICE", &m_EventPhysics);
}


///////////////////////////////////////////////////////////////////////////
void TVOICEPhysics::InitializeRootOutput() {
	TTree* outputTree = RootOutput::getInstance()->GetTree();
	TFile* outputFile = RootOutput::getInstance()->GetFile();
	outputFile->WriteObjectAny(m_EventPhysics,"TVOICEPhysics","VOICE");

/*	outputTree->Branch("nhit",&nhit,"nhit/I");
	outputTree->Branch("type",type,"type[nhit]/I");
	outputTree->Branch("E",E,"E[nhit]/D");
	outputTree->Branch("X",X,"X[nhit]/D");
	outputTree->Branch("Y",Y,"Y[nhit]/D");
	outputTree->Branch("Z",Z,"Z[nhit]/D");
	outputTree->Branch("Time",Time,"Time[nhit]/D");
	outputTree->Branch("Theta",Theta,"Theta[nhit]/D");
	outputTree->Branch("Phi",Phi,"Phi[nhit]/D");
	outputTree->Branch("Atom_Z",Atom_Z,"Atom_Z[nhit]/I");
	outputTree->Branch("A",A,"A[nhit]/I");
	outputTree->Branch("PDG",PDG,"PDG[nhit]/I");*/
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TVOICEPhysics::Construct() {
	return (NPL::VDetector*) new TVOICEPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
	class proxy_VOICE{
		public:
			proxy_VOICE(){
				NPL::DetectorFactory::getInstance()->AddToken("VOICE","VOICE");
				NPL::DetectorFactory::getInstance()->AddDetector("VOICE",TVOICEPhysics::Construct);
			}
	};

	proxy_VOICE p_VOICE;
}

