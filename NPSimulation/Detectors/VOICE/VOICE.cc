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
 *  This class describe  VOICE simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4Colour.hh"
#include "G4Material.hh"

// NPTool header
#include "VOICE.hh"
#include "VOICEScorers.hh"
//#include "CalorimeterScorers.hh"
//#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
#include "NPCore.h"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
using namespace VOICE_NS;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// VOICE Specific Method
VOICE::VOICE(){
	HCID_Anode = -1;
	HCID_Cathode= -1;
	HCID_FG= -1;
	HCID_END= -1;
	HCID_Gas= -1;
	m_Event = new TVOICEData() ;
//	m_EventGas = new TVOICEGasData() ;
	m_AnodeDetector = 0;
	m_FGDetector = 0;
	m_CathodeDetector = 0;
	m_EndDetector = 0;
	m_GasDetector = 0;
	
	m_logicGas=0;
	m_ReactionRegion=NULL;
	m_GasStepSize= 0.5*mm;
	// RGB Color + Transparency
	m_VisSquare 	= new G4VisAttributes(G4Colour(0.8, 0, 0.5));   
	m_VisWire 	= new G4VisAttributes(G4Colour(1, 0, 0));   
	m_VisPCB 	= new G4VisAttributes(G4Colour(0.,0.5, 0.5,0.9));   
	m_VisPCBFG 	= new G4VisAttributes(G4Colour(0.,0.2, 0.8,0.9));   
	m_VisPCBCa 	= new G4VisAttributes(G4Colour(0.8,0.2, 0.2,0.9));   
	m_VisEND 	= new G4VisAttributes(G4Colour(0.5,0.1, 0.5));   
	m_VisGas 	= new G4VisAttributes(G4Colour(0.5,0.2,0.3,0.99));   
	m_cat_zpos.clear();

}

VOICE::~VOICE(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*void VOICE::AddDetector(G4ThreeVector POS, string  Shape){
// Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
m_R.push_back(POS.mag());
m_Theta.push_back(POS.theta());
m_Phi.push_back(POS.phi());
m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void VOICE::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
m_R.push_back(R);
m_Theta.push_back(Theta);
m_Phi.push_back(Phi);
m_Shape.push_back(Shape);
*/

void VOICE::AddDetector(string Type, double zpos, double Tilt, double Rot, double radi, double pitch){
	m_Type.push_back(Type);
	m_zpos.push_back(zpos);
	m_Tilt.push_back(Tilt);
	m_Rot.push_back(Rot);
	m_radi.push_back(radi);
	m_pitch.push_back(pitch);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* VOICE::BuildAnodeDetector(double radi, double pitch){
	if(!m_Anode){
		G4Tubs* Anodewire = new G4Tubs("VOICE_Anode_wire",0,radi,VOICE_NS::WireLength*0.5, 0, 360*deg);
		G4Box* Anodelayer = new G4Box("VOICE_Anodelayer",VOICE_NS::ElectrodePCBX*0.5,VOICE_NS::ElectrodePCBY*0.5,VOICE_NS::ElectrodePCBZ*0.5);
		G4Box* AnodeSub = new G4Box("VOICE_AnodeSub",VOICE_NS::ElectrodeSubX*0.5,VOICE_NS::ElectrodeSubY*0.5,VOICE_NS::ElectrodeSubZ*0.5);
		G4VSolid* AnodeSolidPCB = new G4SubtractionSolid("AnodeSolidPCB",Anodelayer,AnodeSub);

		G4Material* WireMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(VOICE_NS::wireMaterial);
		G4Material* PCBMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(VOICE_NS::PCBMaterial);
		G4LogicalVolume* AnodeWireLogic = new G4LogicalVolume(Anodewire,WireMaterial,"AnodeWireLogic",0,0,0);
		G4LogicalVolume* AnodePCBLogic = new G4LogicalVolume(AnodeSolidPCB,PCBMaterial,"AnodePCBLogic",0,0,0);

		AnodeWireLogic->SetSensitiveDetector(m_AnodeDetector);
		AnodeWireLogic->SetVisAttributes(m_VisWire);
		AnodePCBLogic->SetVisAttributes(m_VisPCB);
		AnodePCBLogic->SetSensitiveDetector(m_AnodeDetector);

		m_Anode = new G4AssemblyVolume();

		G4RotationMatrix* wireRot = new G4RotationMatrix(0,0,0);
		wireRot->rotateX(90*deg);

		G4ThreeVector pos = G4ThreeVector(0,0,0);
		m_Anode->AddPlacedVolume(AnodePCBLogic,pos,0);
		int n_wire = (ceil(VOICE_NS::ElectrodeSubX/pitch));
		double x_start = -n_wire*0.5*pitch;

		    for(int i=0; i<n_wire; i++){
			pos = G4ThreeVector(x_start+i*pitch,0,0);
			m_Anode->AddPlacedVolume(AnodeWireLogic,pos,wireRot); //wirehere
		  }
		//    m_Anode->SetVisAttributes(m_VisSquare);
	}
	return m_Anode;
}

G4AssemblyVolume* VOICE::BuildCathodeDetector(double radi, double pitch){
	if(!m_Cathode){
		G4Tubs* Cathodewire = new G4Tubs("VOICE_Cathode_wire",0,radi,VOICE_NS::WireLength*0.5, 0, 360*deg);
		G4Box* Cathodelayer = new G4Box("VOICE_Cathodelayer",VOICE_NS::ElectrodePCBX*0.5,VOICE_NS::ElectrodePCBY*0.5,VOICE_NS::ElectrodePCBZ*0.5);
		G4Box* CathodeSub = new G4Box("VOICE_CathodeSub",VOICE_NS::ElectrodeSubX*0.5,VOICE_NS::ElectrodeSubY*0.5,VOICE_NS::ElectrodeSubZ*0.5);
		G4VSolid* CathodeSolidPCB = new G4SubtractionSolid("CathodeSolidPCB",Cathodelayer,CathodeSub);

		G4Material* WireMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(VOICE_NS::wireMaterial);
		G4Material* PCBMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(VOICE_NS::PCBMaterial);
		G4LogicalVolume* CathodeWireLogic = new G4LogicalVolume(Cathodewire,WireMaterial,"CathodeWireLogic",0,0,0);
		G4LogicalVolume* CathodePCBLogic = new G4LogicalVolume(CathodeSolidPCB,PCBMaterial,"CathodePCBLogic",0,0,0);

		CathodeWireLogic->SetSensitiveDetector(m_CathodeDetector);
		CathodeWireLogic->SetVisAttributes(m_VisWire);
		CathodePCBLogic->SetVisAttributes(m_VisPCBCa);
		CathodePCBLogic->SetSensitiveDetector(m_CathodeDetector);

		m_Cathode = new G4AssemblyVolume();

		G4ThreeVector pos = G4ThreeVector(0,0,0);
		m_Cathode->AddPlacedVolume(CathodePCBLogic,pos,0);
		int n_wire = (ceil(VOICE_NS::ElectrodeSubX/pitch));
		double x_start = -n_wire*0.5*pitch;

		G4RotationMatrix* wireRot = new G4RotationMatrix(0,0,0);
		wireRot->rotateX(90*deg);
		for(int i=0; i<n_wire; i++){
			pos = G4ThreeVector(x_start+i*pitch,0,0);
			m_Cathode->AddPlacedVolume(CathodeWireLogic,pos,wireRot); //wirehere
		}
		//    m_Cathode->SetVisAttributes(m_VisSquare);
	}
	return m_Cathode;
}

G4AssemblyVolume* VOICE::BuildFGDetector(double radi, double pitch){
	if(!m_FG){
		G4Tubs* FGwire = new G4Tubs("VOICE_FG_wire",0,radi,VOICE_NS::WireLength*0.5, 0, 360*deg);
		G4Box* FGlayer = new G4Box("VOICE_FGlayer",VOICE_NS::ElectrodePCBX*0.5,VOICE_NS::ElectrodePCBY*0.5,VOICE_NS::ElectrodePCBZ*0.5);
		G4Box* FGSub = new G4Box("VOICE_FGSub",VOICE_NS::ElectrodeSubX*0.5,VOICE_NS::ElectrodeSubY*0.5,VOICE_NS::ElectrodeSubZ*0.5);
		G4VSolid* FGSolidPCB = new G4SubtractionSolid("FGSolidPCB",FGlayer,FGSub);

		G4Material* WireMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(VOICE_NS::wireMaterial);
		G4Material* PCBMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(VOICE_NS::PCBMaterial);
		G4LogicalVolume* FGWireLogic = new G4LogicalVolume(FGwire,WireMaterial,"FGWireLogic",0,0,0);
		G4LogicalVolume* FGPCBLogic = new G4LogicalVolume(FGSolidPCB,PCBMaterial,"FGPCBLogic",0,0,0);

		FGWireLogic->SetSensitiveDetector(m_FGDetector);
		FGWireLogic->SetVisAttributes(m_VisWire);
		FGPCBLogic->SetVisAttributes(m_VisPCBFG);
		FGPCBLogic->SetSensitiveDetector(m_FGDetector);

		m_FG = new G4AssemblyVolume();
		G4RotationMatrix* wireRot = new G4RotationMatrix(0,0,0);
		wireRot->rotateX(90*deg);

		G4ThreeVector pos = G4ThreeVector(0,0,0);
		m_FG->AddPlacedVolume(FGPCBLogic,pos,0);
		int n_wire = (ceil(VOICE_NS::ElectrodeSubX/pitch));
		double x_start = -n_wire*0.5*pitch;

		for(int i=0; i<n_wire; i++){
			pos = G4ThreeVector(x_start+i*pitch,0,0);
			m_FG->AddPlacedVolume(FGWireLogic,pos,wireRot);//wirehere
		}
		//    m_FG->SetVisAttributes(m_VisSquare);
	}
	return m_FG;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* VOICE::BuildEndDetector(){
	if(!m_End){

		G4Box* enddet = new G4Box("VOICE_End",VOICE_NS::ElectrodePCBX*0.5,VOICE_NS::ElectrodePCBY*0.5,VOICE_NS::endthick*0.5);
		G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(VOICE_NS::wireMaterial);
		G4LogicalVolume* ENDLogic = new G4LogicalVolume(enddet,DetectorMaterial,"logic_VOICE_End",0,0,0);
		ENDLogic->SetVisAttributes(m_VisEND);
		ENDLogic->SetSensitiveDetector(m_EndDetector);

		m_End = new G4AssemblyVolume();
		G4ThreeVector pos = G4ThreeVector(0,0,0);
		m_End->AddPlacedVolume(ENDLogic,pos,0);

	}
	return m_End;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void VOICE::ReadConfiguration(NPL::InputParser parser){


	vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("VOICE");
	if(NPOptionManager::getInstance()->GetVerboseLevel())
		cout << "//// " << blocks.size() << " Components found"<< endl; 

	//  vector<string> cart = {"POS","Shape"};
	//  vector<string> sphe = {"R","Theta","Phi","Shape"};
	vector<string> voice_vec = {"Type"};
	vector<string> voice_gas_vec = {"gas","pressure"};

	bool FG_flag = false;
	bool Cathode_flag = false;
	for(unsigned int i = 0 ; i < blocks.size() ; i++){
		int count_gas_info =0;
		if(blocks[i]->HasTokenList(voice_vec)){
			if(NPOptionManager::getInstance()->GetVerboseLevel())
				cout<< endl <<"\\\\ VOICE " << i+1 <<endl;

			string Type = blocks[i]->GetString("Type");
			double zpos = -9999;
			if(blocks[i]->HasToken("z-pos"))zpos= blocks[i]->GetDouble("z-pos","mm");
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
			if(strcmp(m_gas.c_str(),"multigas")==0){
				m_pressure = blocks[i]->GetDouble("pressure","atmosphere");
				mulgas_density = blocks[i]->GetDouble("density","g/cm3");
				n_gas = blocks[i]->GetInt("n_gas");
				for(int j=0; j<n_gas; j++){
					string gasstring="gas";
					string fracstring="massfrac";
					gasstring += to_string(j+1);
					fracstring += to_string(j+1);
					m_mulgas[j]=blocks[i]->GetString(gasstring);
					mass_fraction[j]=blocks[i]->GetDouble(fracstring,"");
					cout<<"READ Gas: "<<m_mulgas[j]<<endl;
				}
				count_gas_info++;
			}else{
				//m_pressure = blocks[i]->GetDouble("pressure","torr");
				m_pressure = blocks[i]->GetDouble("pressure","atmosphere");
				count_gas_info++;
			}
		}
		else{
			cout << "ERROR: check your input file formatting " << endl;
			exit(1);
		}
	}
	
	for (unsigned short i = 0 ; i < m_Type.size() ; i++) {
		if(m_Type[i]=="Anode"){
			if(maxZ < m_zpos[i]) maxZ = m_zpos[i];
		}
	}

	total_length = (maxZ+(m_AC_gap*2)+VOICE_NS::ElectrodePCBZ);
	zoffset=(total_length-VOICE_NS::ElectrodePCBZ)*0.5-m_AC_gap;
	
	m_cat_zpos.push_back(m_zpos[0]-m_AC_gap-zoffset);
	for (unsigned short i = 0 ; i < m_Type.size() ; i++) {
		if(m_Type[i]=="Anode"){
			m_cat_zpos.push_back(m_zpos[i]+m_AC_gap-zoffset);
		}
	}
}

G4Material* VOICE::SetGasType(string gas, double pressure){

	G4Material* gasmat = nullptr;
	G4int z=0;
	double filldensity = 0.0;
	int kk=0;
	G4Element* C = new G4Element("Carbon", "C", z=6, 12.011*g/mole);
	G4Element* F = new G4Element("Flourine", "F", z=9, 18.998*g/mole);
	G4Material* m_Ar = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_Ar");
	G4Material* m_He = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_He");
	G4Material* m_H = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_H");
	G4Material* m_methane = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_METHANE");

	if(strcmp(gas.c_str(),"multigas")==0){
		if(n_gas<2 || n_gas>10){
			cout<<"Number of gas components is out of the range"<<endl;
			exit(1);
		}else{
			gasmat=new G4Material("gasmat",(pressure/atmosphere)*(mulgas_density/(g/cm3))*g/cm3,n_gas,kStateGas,300.*kelvin,pressure);
			for(int i=0; i<n_gas; i++){
				G4Material* temp_mat = MaterialManager::getInstance()->GetMaterialFromLibrary(m_mulgas[i].c_str());
				if(temp_mat==NULL){
					cout<<gas<<" is not defined material for Gas type."<<endl;
					exit(1);
				}
				gasmat->AddMaterial(temp_mat,mass_fraction[i]);	
			}
		}
	}else if(strcmp(gas.c_str(),"CF4")==0){
		filldensity = 0.00366*g/cm3;
		gasmat=new G4Material("gasmat",(pressure/atmosphere)*(filldensity/(g/cm3))*g/cm3,2,kStateGas,300.*kelvin,pressure);
		gasmat->AddElement(C,1);
		gasmat->AddElement(F,4);
	}else if(strcmp(gas.c_str(),"P10")==0){
		filldensity = 0.00156*g/cm3;	
		gasmat=new G4Material("gasmat",(pressure/atmosphere)*(filldensity/(g/cm3))*g/cm3,2,kStateGas,300.*kelvin,pressure);
		gasmat->AddMaterial(m_Ar,0.9574);		//mass fraction 
		gasmat->AddMaterial(m_methane,0.0426);
	}else if(strcmp(gas.c_str(),"C3F8")==0){
		filldensity = 0.00817*g/cm3;
		gasmat=new G4Material("gasmat",(pressure/atmosphere)*(filldensity/(g/cm3))*g/cm3,2,kStateGas,300.*kelvin,pressure);
		gasmat->AddElement(C,3);
		gasmat->AddElement(F,8);

	}else{	
		gasmat = MaterialManager::getInstance()->GetMaterialFromLibrary(gas.c_str());
		filldensity = (pressure/atmosphere)*(gasmat->GetDensity()/(g/cm3))*g/cm3;
		gasmat = MaterialManager::getInstance()->GetMaterialFromLibrary(gas.c_str(),filldensity);

		if(gasmat==NULL){
			cout<<gas<<" is not defined material for Gas type."<<endl;
			exit(1);
		}
	}

	/*
	   else if(strcmp(gas.c_str(),"CF4")==0){
	   filldensity = 0.00366*g/cm3;
	   gasmat=new G4Material("gasmat",(pressure/atmosphere)*(filldensity/(g/cm3))*g/cm3,2,kStateGas,300.*kelvin,pressure);
	   gasmat->AddElement(C,1);
	   gasmat->AddElement(F,4);
	   }else if(strcmp(gas.c_str(),"P10")==0){
	   filldensity = 0.00156*g/cm3;
	   gasmat=new G4Material("gasmat",(pressure/atmosphere)*(filldensity/(g/cm3))*g/cm3,2,kStateGas,300.*kelvin,pressure);
	   gasmat->AddMaterial(m_Ar,0.9574);
	   gasmat->AddMaterial(m_methane,0.0426);
	   }else if(strcmp(gas.c_str(),"C3F8")==0){
	   filldensity = 0.00817*g/cm3;
	   gasmat=new G4Material("gasmat",(pressure/atmosphere)*(filldensity/(g/cm3))*g/cm3,2,kStateGas,300.*kelvin,pressure);
	   gasmat->AddElement(C,3);
	   gasmat->AddElement(F,8);
	   }else if(strcmp(gas.c_str(),"He")==0){
	   filldensity = 0.00017*g/cm3;
	   gasmat=new G4Material("gasmat",(pressure/atmosphere)*(filldensity/(g/cm3))*g/cm3,1,kStateGas,300.*kelvin,pressure);
	   gasmat->AddMaterial(m_He,1.0);
	   }else{
	//		gasmat = G4Material::GetMaterial(gas);
	gasmat = MaterialManager::getInstance()->GetMaterialFromLibrary(gas);
	filldensity = (pressure/atmosphere)*(gasmat->GetDensity()/(g/cm3))*g/cm3;
	gasmat = MaterialManager::getInstance()->GetMaterialFromLibrary(gas,filldensity);

	if(gasmat==NULL){
	cout<<gas<<" is not defined material for Gas type."<<endl;
	}
	}*/

	cout<<gasmat->GetChemicalFormula()<<"\t"<<gasmat->GetPressure()/atmosphere<<"\t"<<gasmat->GetDensity()/(g/cm3)<<" is the Gas."<<endl;
	return gasmat;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void VOICE::ConstructDetector(G4LogicalVolume* world){

	G4Material* VOICE_gas = SetGasType(m_gas,m_pressure);

	if(VOICE_gas==NULL) {
		cout<<"ERROR: Wrong gas input. Exit."<<endl;
		exit(1);
	}

	//	world->SetMaterial(VOICE_gas);

	int countdet=0;
	vector<double> anode_zpos;

	// edit by minju
	G4Box* GasBox = new G4Box("VOICE_GasBox",VOICE_NS::ElectrodePCBX*0.5,VOICE_NS::ElectrodePCBY*0.5,total_length*0.5+VOICE_NS::GasDeadLayer);
	//double zoffset=(maxZ+VOICE_NS::endthick+VOICE_NS::ElectrodeSubZ)*0.5;
	m_logicGas = new G4LogicalVolume(GasBox,VOICE_gas,"logicTarget");
	m_logicGas->SetVisAttributes(m_VisGas);
	m_logicGas->SetSensitiveDetector(m_GasDetector);
	//	m_logicGas->SetVisAttributes(G4VisAttributes::GetInvisible());
	if(!m_ReactionRegion){
		G4ProductionCuts* ecut = new G4ProductionCuts();
		ecut->SetProductionCut(1000*mm,"e-");
		m_ReactionRegion= new G4Region("NPSimulationProcess");
		m_ReactionRegion->SetProductionCuts(ecut);
		m_ReactionRegion->AddRootLogicalVolume(m_logicGas);
		//		m_ReactionRegion->AddRootLogicalVolume(m_logicGas);
		//		m_ReactionRegion->AddRootLogicalVolume(m_logicGas);
		//		cout<<"Number of ROOTVOLUMES: "<<m_ReactionRegion->GetNumberOfRootVolumes()<<endl;
		m_ReactionRegion->SetUserLimits(new G4UserLimits(m_GasStepSize));
	}
	G4FastSimulationManager* mng = m_ReactionRegion->GetFastSimulationManager();
	unsigned int size = m_ReactionModel.size();
	for(unsigned int i = 0 ; i < size ; i++){
		mng->RemoveFastSimulationModel(m_ReactionModel[i]);
	}
	m_ReactionModel.clear();
	G4VFastSimulationModel* fsm;
	fsm = new NPS::BeamReaction("BeamReaction",m_ReactionRegion);
	m_ReactionModel.push_back(fsm);
	fsm = new NPS::Decay("Decay",m_ReactionRegion);
	m_ReactionModel.push_back(fsm);

	//	G4Region* Region_cut = new G4Region("RegionCut");
	//	Region_cut->SetProductionCuts(ecut);
	//	Region_cut->AddRootLogicalVolume(m_SquareDetector);

	//G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0),m_logicGas,"VOICE_Gas",world,false,0);

	for (unsigned short i = 0 ; i < m_Type.size() ; i++) {

		// So the face of the detector is at R instead of the middle
		// Det_pos+=Det_pos.unit()*VOICE_NS::Thickness*0.5;

		cout<<"========== STARTING DETECTOR CONSTRUCTION ==========="<<endl;
		G4double wX = 0.0;
		G4double wY = 0.0;
		G4double wZ = 0.0;   
		G4AssemblyVolume*det;

		if(m_Type[i]=="Anode"){
			det = BuildAnodeDetector(m_radi[i],m_pitch[i]);
			anode_zpos.push_back(m_zpos[i]);
			cout<<"========== ANODE IS FABRICATED ==========="<<endl;
		}
		else if(m_Type[i]=="FG") {
			det = BuildFGDetector(m_radi[i],m_pitch[i]);

			cout<<"========== FG IS FABRICATED ==========="<<endl;}
		else if(m_Type[i]=="Cathode"){ det = BuildCathodeDetector(m_radi[i],m_pitch[i]);
			cout<<"========== CATHODE IS FABRICATED ==========="<<endl;}
		else if(m_Type[i]=="END"){ det = BuildEndDetector();
			cout<<"========== END DETECTOR IS FABRICATED ==========="<<endl;}
		else{
			std::cerr<<"no type "<<m_Type[i]<< " exists."<<std::endl;
			continue;	
		}

		G4ThreeVector Det_pos = G4ThreeVector(0, 0, 0) ;
		G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
		G4double wTilt = m_Tilt[i];
		G4double wRot = m_Rot[i];
		Rot->rotateZ(wRot*deg);
		Rot->rotateX(wTilt*deg);

		if(m_Type[i]=="Anode"){
			wZ = m_zpos[i];
			Det_pos=G4ThreeVector(wX,wY,wZ-zoffset);
			det->MakeImprint(m_logicGas,Det_pos,Rot,countdet+1,false);
			//det->MakeImprint(world,Det_pos,Rot,countdet+1);
			countdet++;
			cout<<"========== ANODE"<<countdet<<" IS PLACED ==========="<<endl;
		}else if(m_Type[i]=="FG"){
			for(int j=0; j<anode_zpos.size(); j++){
				Det_pos = G4ThreeVector(wX,wY,anode_zpos[j]-m_AFG_gap-zoffset);
				det->MakeImprint(m_logicGas,Det_pos,Rot,countdet+1,false);
				//det->MakeImprint(world,Det_pos,Rot,countdet+1);
				countdet++;
				Det_pos = G4ThreeVector(wX,wY,anode_zpos[j]+m_AFG_gap-zoffset);
				det->MakeImprint(m_logicGas,Det_pos,Rot,countdet+1,false);
				//det->MakeImprint(world,Det_pos,Rot,countdet+1);
				countdet++;
			}
			cout<<"========== FG IS PLACED ==========="<<endl;
		}else if(m_Type[i]=="Cathode"){
			
			Det_pos = G4ThreeVector(wX,wY,anode_zpos[0]-m_AC_gap-zoffset);
			det->MakeImprint(m_logicGas,Det_pos,Rot,countdet+1,false);		
			countdet++;
			
			for(int j=0; j<anode_zpos.size(); j++){
				Det_pos = G4ThreeVector(wX,wY,anode_zpos[j]+m_AC_gap-zoffset);
				det->MakeImprint(m_logicGas,Det_pos,Rot,countdet+1,false);		
				//det->MakeImprint(world,Det_pos,Rot,countdet+1);
				countdet++;
			}
			cout<<"========== CATHODE IS PLACED ==========="<<endl;
		}else if(m_Type[i]=="END"){
			if(i<m_Type.size()-1){
				cout<<"@@@@@@@@@ ERROR: Put END detector information at the last"<<endl;
				exit(1);
			}
			cout<<"maxZ"<<maxZ<<endl;
			//edit by minju
			wZ = maxZ+m_AC_gap+endthick*0.5+VOICE_NS::GasDeadLayer;
			//			wZ = maxZ+endthick*0.5+m_AFG_gap;
			Det_pos = G4ThreeVector(wX,wY,wZ-zoffset);
			det->MakeImprint(world,Det_pos,Rot,countdet+1,false);
			//det->MakeImprint(world,Det_pos,Rot,countdet+1);
			countdet++;
			cout<<"========== END DETECTOR IS PLACED ==========="<<endl;
		}

		//iterator is equal to fPVStore.begin()
		std::vector<G4VPhysicalVolume*>::iterator it = det ->GetVolumesIterator();
		unsigned int NbrImprints = det->GetImprintsCount();
		unsigned int NbrTotalPV = det->TotalImprintedVolumes();
		unsigned int NbrComponents = NbrTotalPV/NbrImprints;
		// set copy numbers of components of assembly volume to the current detector number
		for(it += (NbrImprints-1)*NbrComponents; it<= det->GetVolumesIterator()+NbrTotalPV-1; it++)
			(*it)->SetCopyNo(i+1);
	}

	//	new G4PVPlacement(0,G4ThreeVector(0,0,0),m_logicStage,"VOICE_Stage_Volume",world,true,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,0),m_logicGas,"VOICE_Gas_Volume",world,true,0);

	cout<<"========== TOTAL "<<countdet<<" COMPONENTS ARE PLACED ==========="<<endl;
	std::cout<<"Detector Construction is completed"<<std::endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void VOICE::InitializeRootOutput(){
	RootOutput *pAnalysis = RootOutput::getInstance();
	TTree *pTree = pAnalysis->GetTree();
	if(!pTree->FindBranch("VOICE")){
		pTree->Branch("VOICE", "TVOICEData", &m_Event) ;
	}
	pTree->SetBranchAddress("VOICE", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void VOICE::ReadSensitive(const G4Event* event){
	m_Event->Clear();
	//	m_EventGas->Clear();
	//////////////////May start modification here!!!! Sunghan!!!
	auto HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	if(HCID_Anode == -1)
		HCID_Anode = G4SDManager::GetSDMpointer()->GetCollectionID("AnodeDet/AnodeScorer");
	if(HCID_Cathode == -1)
		HCID_Cathode = G4SDManager::GetSDMpointer()->GetCollectionID("CathodeDet/CathodeScorer");
	if(HCID_FG == -1)
		HCID_FG = G4SDManager::GetSDMpointer()->GetCollectionID("FGDet/FGScorer");
	if(HCID_END == -1)
		HCID_END = G4SDManager::GetSDMpointer()->GetCollectionID("ENDDet/ENDScorer");
	if(HCID_Gas == -1)
		HCID_Gas = G4SDManager::GetSDMpointer()->GetCollectionID("GasDet/GasScorer");


	// loop for the event map

	auto evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_Gas));
	map<G4int, G4double**>::iterator it;
//	cout<<evtMap->GetSize()<<"asdfasdfasdfasdfasdf"<<endl;
	for(it = evtMap-> GetMap()->begin(); it!= evtMap->GetMap()->end(); it++){
		m_Event->Set(4,
				(*(it->second))[0], //init kin. energy
				(*(it->second))[1], //init x
				(*(it->second))[2], //init y
				(*(it->second))[3], //init z
				(*(it->second))[4], //time
				(*(it->second))[5], //theta
				(*(it->second))[6], //phi
				(*(it->second))[7], //atomic number (Z)
				(*(it->second))[8],   //atomic mass (A)
				(*(it->second))[9],   //PDG encoding
				(*(it->second))[10], //current kin.energy
				(*(it->second))[11], //xf
				(*(it->second))[12], //yf
				(*(it->second))[13], //zf
				(*(it->second))[14], //timef
				(*(it->second))[15], //thetaf
				(*(it->second))[16] //phif
			    );
	}

	evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_Anode));
	for(it = evtMap-> GetMap()->begin(); it!= evtMap->GetMap()->end(); it++){
		m_Event->Set(0,
				(*(it->second))[0], //energy deposition
				(*(it->second))[1], //hit x
				(*(it->second))[2], //hit y
				(*(it->second))[3], //hit z
				(*(it->second))[4], //time
				(*(it->second))[5], //theta
				(*(it->second))[6], //phi
				(*(it->second))[7], //atomic number (Z)
				(*(it->second))[8],   //atomic mass (A)
				(*(it->second))[9],   //PDG encoding
				-1,-1,-1,-1,-1,-1,-1
			    );
	}

	evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_Cathode));
	for(it = evtMap-> GetMap()->begin(); it!= evtMap->GetMap()->end(); it++){
		m_Event->Set(1,
				(*(it->second))[0], //energy deposition
				(*(it->second))[1], //hit x
				(*(it->second))[2], //hit y
				(*(it->second))[3], //hit z
				(*(it->second))[4], //time
				(*(it->second))[5], //theta
				(*(it->second))[6], //phi
				(*(it->second))[7], //atomic number (Z)
				(*(it->second))[8],   //atomic mass (A)
				(*(it->second))[9],   //PDG encoding
				-1,-1,-1,-1,-1,-1,-1
			    );
	}

	evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_FG));
	for(it = evtMap-> GetMap()->begin(); it!= evtMap->GetMap()->end(); it++){
		m_Event->Set(2,
				(*(it->second))[0], //energy deposition
				(*(it->second))[1], //hit x
				(*(it->second))[2], //hit y
				(*(it->second))[3], //hit z
				(*(it->second))[4], //time
				(*(it->second))[5], //theta
				(*(it->second))[6], //phi
				(*(it->second))[7], //atomic number (Z)
				(*(it->second))[8],   //atomic mass (A)
				(*(it->second))[9],   //PDG encoding
				-1,-1,-1,-1,-1,-1,-1
			    );
	}

	evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_END));
	for(it = evtMap-> GetMap()->begin(); it!= evtMap->GetMap()->end(); it++){
		m_Event->Set(3,
				(*(it->second))[0], //energy deposition
				(*(it->second))[1], //hit x
				(*(it->second))[2], //hit y
				(*(it->second))[3], //hit z
				(*(it->second))[4], //time
				(*(it->second))[5], //theta
				(*(it->second))[6], //phi
				(*(it->second))[7], //atomic number (Z)
				(*(it->second))[8],   //atomic mass (A)
				(*(it->second))[9],   //PDG encoding
				-1,-1,-1,-1,-1,-1,-1
			    );
	}



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void VOICE::InitializeScorers() { 
	// This check is necessary in case the geometry is reloaded
	bool already_exist_A = false; 
	bool already_exist_C = false; 
	bool already_exist_F = false; 
	bool already_exist_E = false; 
	bool already_exist_G = false; 
	m_AnodeDetector = CheckScorer("AnodeDet",already_exist_A) ;
	m_CathodeDetector = CheckScorer("CathodeDet",already_exist_C) ;
	m_FGDetector = CheckScorer("FGDet",already_exist_F) ;
	m_EndDetector = CheckScorer("ENDDet",already_exist_E) ;
	m_GasDetector = CheckScorer("GasDet",already_exist_G) ;
	// Otherwise the scorer is initialised
	// See NPSimultion/Scorers/VOICEScorer.xx
	if(!already_exist_A){
		G4VPrimitiveScorer* AnodeScorer = new VOICEScorers::PS_Electrode("AnodeScorer",0,ElectrodeSubX,ElectrodeSubY,ElectrodePCBZ,0);
		m_AnodeDetector->RegisterPrimitive(AnodeScorer);
		G4SDManager::GetSDMpointer()->AddNewDetector(m_AnodeDetector) ;
	}

	if(!already_exist_C){
		G4VPrimitiveScorer* CathodeScorer = new VOICEScorers::PS_Electrode("CathodeScorer",0,ElectrodeSubX,ElectrodeSubY,ElectrodePCBZ,0);
		m_CathodeDetector->RegisterPrimitive(CathodeScorer);
		G4SDManager::GetSDMpointer()->AddNewDetector(m_CathodeDetector) ;
	}

	if(!already_exist_F){
		G4VPrimitiveScorer* FGScorer = new VOICEScorers::PS_Electrode("FGScorer",0,ElectrodeSubX,ElectrodeSubY,ElectrodePCBZ,0);
		m_FGDetector->RegisterPrimitive(FGScorer);
		G4SDManager::GetSDMpointer()->AddNewDetector(m_FGDetector) ;
	}

	if(!already_exist_E){
		G4VPrimitiveScorer* ENDScorer = new VOICEScorers::PS_END("ENDScorer",0,endthick,0);
		m_EndDetector->RegisterPrimitive(ENDScorer);
		G4SDManager::GetSDMpointer()->AddNewDetector(m_EndDetector) ;
	}
	if(!already_exist_G){
		G4VPrimitiveScorer* GasScorer = new VOICEScorers::PS_Gas("GasScorer",m_cat_zpos,m_GasStepSize,0,0);
		m_GasDetector->RegisterPrimitive(GasScorer);
		G4SDManager::GetSDMpointer()->AddNewDetector(m_GasDetector) ;
	}


	cout<<"//////////////////////// INITIALIZING SCORER IS DONE //////////////////////////////"<<endl<<endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* VOICE::Construct(){
	return  (NPS::VDetector*) new VOICE();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
	class proxy_nps_VOICE{
		public:
			proxy_nps_VOICE(){
				NPS::DetectorFactory::getInstance()->AddToken("VOICE","VOICE");
				NPS::DetectorFactory::getInstance()->AddDetector("VOICE",VOICE::Construct);
			}
	};

	proxy_nps_VOICE p_nps_VOICE;
}
