#ifndef NPBreitWigner_h
#define NPBreitWigner_h
#if __cplusplus > 201703L
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author : Audrey ANNE  contact address: anne@lpccaen.in2p3.fr     *
 *                                                                           *
 * Creation Date   : February 2024                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>

// ROOT header
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include <cmath>
#include "TRandom.h"

// NPL header
#include "NPInputParser.h"
#include "NPParticle.h"

class TF1;

using namespace std;

namespace NPL{
  
  class BreitWigner{
    
  public:  // Constructors and Destructors
    BreitWigner();
    ~BreitWigner();
 

  public: //Read Method
    void ReadConfigurationFile(string Path);
    void ReadConfigurationFile(NPL::InputParser);
  

  private:
    //Physical Constants
    const double hc;
    const double uma;
    
    bool fShiftZero; // = 0 -> Shift ; = 1 -> no Shift 
    bool fWidthCst; // = 0 -> constant Width ; = 1 -> Width 

    double fR0;

    double fX0;
    double fW0;
    int fell; // l
    int fA;  // Fragment's A number
    int fn;  // Neutrons number
    double fparam0;

    Particle fHeavy; // Heavy in BreitWigner (example : He8 -> He4 + 4n, Heavy is He8).
    Particle fLight;  // Light in BreitWigner (example : He8 -> He4 + 4n, Light is He4).
    Particle fNeutron;

    //Calculated values
    double frho_0;
    double fP0;
    double fSF0;
    double fmu; //
    double fR;

    double fQ;

 
  public:

    void SetA(int A){fA = A;}
    int GetA()const{return fA;}

    void SetL(int l){fell = l;}
    int GetL()const{return fell;}

    double Calculate_Rho_Var(double E) const;
    double Width(double rho_X) const;
    double Shift(double rho_X) const;
    double Width_Ratio(double rho_X) const;
    double Shift_Ratio(double rho_x) const;
    double Shift_Factor(double rho) const;
    double Penetration_Factor(double rho) const;

    
    double Eval(const double & x_) const;
    double Eval_For_TF1(const double * x_, const double * /* params_ */) const;

    
    TF1 * Make_TF1() const;
    TH1D * Make_TH1() const;

    //TF1 * Make_TF1(const char * name_) const;
    
    //only for test -> TF1 and TH1 put in .root file
    void Save_Histo() const ;
    void Save_TF1() const ;
    
    void Dump() const;
  }; 
}

#endif //c++17
#endif //ndef
