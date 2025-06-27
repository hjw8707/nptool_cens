#ifndef __VOICEDATA__
#define __VOICEDATA__
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

// STL
#include <vector>
using namespace std;

// ROOT
#include "TObject.h"

class TVOICEData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // Energy
	vector<Int_t>		fType;
	vector<Double_t>	fE;
	vector<Double_t>	fX;
	vector<Double_t>	fY;
	vector<Double_t>	fZ;
	vector<Double_t>	fTime;
	vector<Double_t>	fTheta;
	vector<Double_t>	fPhi;
	vector<Int_t>		fAtom_Z;
	vector<Int_t>		fA;
	vector<Int_t>		fPDG;
	vector<Double_t>	fE_f;
	vector<Double_t>	fX_f;
	vector<Double_t>	fY_f;
	vector<Double_t>	fZ_f;
	vector<Double_t>	fTime_f;
	vector<Double_t>	fTheta_f;
	vector<Double_t>	fPhi_f;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TVOICEData();
    ~TVOICEData();
    

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public:
    void Clear();
    void Clear(const Option_t*) {};
    void Dump() const;


  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////
    // Energy
    inline void Set(const Int_t Type, const Double_t Energy, const Double_t X, const Double_t Y, const Double_t Z,
		    const Double_t Time, const Double_t Theta, const Double_t Phi, const Int_t Atom_Z, const Int_t A, const Int_t PDG, 
		    const Double_t Energy_f, const Double_t X_f, const Double_t Y_f,const Double_t Z_f,const Double_t Time_f, const Double_t Theta_f, const Double_t Phi_f){
	fType.push_back(Type);
	fE.push_back(Energy);
	fX.push_back(X);
	fY.push_back(Y);
	fZ.push_back(Z);
	fTime.push_back(Time);
	fTheta.push_back(Theta);
	fPhi.push_back(Phi);
	fAtom_Z.push_back(Atom_Z);
	fA.push_back(A);
	fPDG.push_back(PDG);
	fE_f.push_back(Energy_f);
	fX_f.push_back(X_f);
	fY_f.push_back(Y_f);
	fZ_f.push_back(Z_f);
	fTime_f.push_back(Time_f);
	fTheta_f.push_back(Theta_f);
	fPhi_f.push_back(Phi_f);
    }

    //////////////////////    GETTERS    ////////////////////////
	inline Int_t GetMult() const {return fType.size();}
	inline Int_t GetType(Int_t i) const {return fType[i];}
	inline Double_t GetEnergy(Int_t i) const {return fE[i];}
	inline Double_t GetX(Int_t i) const {return fX[i];}
	inline Double_t GetY(Int_t i) const {return fY[i];}
	inline Double_t GetZ(Int_t i) const {return fZ[i];}
	inline Double_t GetTime(Int_t i) const {return fTime[i];}
	inline Double_t GetTheta(Int_t i) const {return fTheta[i];}
	inline Double_t GetPhi(Int_t i) const {return fPhi[i];}
	inline Int_t GetAtom_Z(Int_t i) const {return fAtom_Z[i];}
	inline Int_t GetA(Int_t i) const {return fA[i];}
	inline Int_t GetPDG(Int_t i) const {return fPDG[i];}
	inline Double_t GetEnergy_f(Int_t i) const {return fE_f[i];}
	inline Double_t GetX_f(Int_t i) const {return fX_f[i];}
	inline Double_t GetY_f(Int_t i) const {return fY_f[i];}
	inline Double_t GetZ_f(Int_t i) const {return fZ_f[i];}
	inline Double_t GetTime_f(Int_t i) const {return fTime_f[i];}
	inline Double_t GetTheta_f(Int_t i) const {return fTheta_f[i];}
	inline Double_t GetPhi_f(Int_t i) const {return fPhi_f[i];}

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TVOICEData,1)  // VOICEData structure
};

#endif
