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
 *  This class hold VOICE Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TVOICEData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TVOICEData)

//////////////////////////////////////////////////////////////////////
TVOICEData::TVOICEData() {
}

//////////////////////////////////////////////////////////////////////
TVOICEData::~TVOICEData() {
}

//////////////////////////////////////////////////////////////////////
void TVOICEData::Clear() {
	fType.clear();
	fE.clear();
	fX.clear();
	fY.clear();
	fZ.clear();
	fTime.clear();

	fTheta.clear();
	fPhi.clear();
	fAtom_Z.clear();
	fA.clear();
	fPDG.clear();
	
	fE_f.clear();
	fX_f.clear();
	fY_f.clear();
	fZ_f.clear();
	fTheta_f.clear();
	fPhi_f.clear();
	fTime_f.clear();
}



//////////////////////////////////////////////////////////////////////
void TVOICEData::Dump() const {
  // This method is very useful for debuging and worth the dev.
	std::cout << "=========== Check VOICE Data ==========" << std::endl;
	std::cout << " Total size = " <<GetMult() << std::endl;

	for(Int_t i=0; i<GetMult(); i++){
		std::cout<<" Type = "<<fType[i]<<", ";
		std::cout<<" Energy = "<<fE[i]<<", ";
		std::cout<<" x-pos = "<<fX[i]<<", ";
		std::cout<<" y-pos = "<<fY[i]<<", ";
		std::cout<<" z-pos = "<<fZ[i]<<", ";
		std::cout<<" Time = "<<fTime[i]<<", ";
		std::cout<<" Theta = "<<fTheta[i]<<", ";
		std::cout<<" Phi = "<<fPhi[i]<<", ";
		std::cout<<" At. Num = "<<fAtom_Z[i]<<", ";
		std::cout<<" Mass = "<<fA[i]<<std::endl;
		std::cout<<" PDG = "<<fPDG[i]<<std::endl;
		std::cout<<" Energy_f = "<<fE_f[i]<<", ";
		std::cout<<" x-pos_f = "<<fX_f[i]<<", ";
		std::cout<<" y-pos_f = "<<fY_f[i]<<", ";
		std::cout<<" z-pos_f = "<<fZ_f[i]<<", ";
		std::cout<<" Time_f = "<<fTime_f[i]<<", ";
		std::cout<<" Theta_f = "<<fTheta_f[i]<<", ";
		std::cout<<" Phi_f = "<<fPhi_f[i]<<", ";
	}
	std::cout << "====================================" << std::endl;
}
