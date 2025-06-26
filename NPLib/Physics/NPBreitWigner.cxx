#if __cplusplus > 201703L
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Audrey ANNE contact address: anne@lpccaen.in2p3.fr     *
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

// STL
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <vector>

// NPL header
#include "NPBreitWigner.h"
#include "NPFunction.h"
#include "NPOptionManager.h"

// ROOT Header
#include "TDirectory.h"
#include "TCanvas.h"

// Use CLHEP System of unit and Physical Constant
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"

#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

#include "TAxis.h"

using namespace NPL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
BreitWigner::BreitWigner():hc(hbarc/fermi), uma(amu_c2){
  
  fR0 = 1.4; // Value in Breit Wigner models
  fNeutron = Particle("n");

  fShiftZero = 0;
  fWidthCst = 0;

  // Breit Wigner parameters
  fX0 = 1.0;
  fW0 = 1.0;
  fell = 0;
  fA = 1;
  fn = 1;
  fparam0 = 1.0;
  
  // Calculated values
  fmu = 1.0;
  fR = 1.0;
  frho_0 = 1.0;
  fP0 = 1.0;
  fSF0 = 1.0;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BreitWigner::~BreitWigner(){
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BreitWigner::ReadConfigurationFile(string Path){
  NPL::InputParser parser(Path);
  ReadConfigurationFile(parser);
}

void BreitWigner::ReadConfigurationFile(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("BreitWigner");
  if(blocks.size()>0 && NPOptionManager::getInstance()->GetVerboseLevel())
    cout << endl << "\033[1;35m//// Breit Wigner energy distribution " << endl; 

  vector<string> token   = {"Energy0","Width0", "l", "Heavy","Light", "NeutronNumber","MultiplicateurAmplitude", "ShiftZero", "WidthCst"};
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
 
      fX0 = blocks[i]->GetDouble("Energy0","MeV");
      fW0 = blocks[i]->GetDouble("Width0","MeV");
      if(fW0<=0) fW0 = 1e-4 ;
      fell = blocks[i]->GetInt("l");
      fHeavy = Particle(blocks[i]->GetString("Heavy"));
      fLight =  Particle(blocks[i]->GetString("Light"));
      fn = blocks[i]->GetInt("NeutronNumber");
      fparam0 =  blocks[i]->GetDouble("MultiplicateurAmplitude","MeV");
      fShiftZero =  blocks[i]->GetInt("ShiftZero");
      fWidthCst =  blocks[i]->GetInt("WidthCst");

    }

    else{
      cout << "ERROR: check your input file formatting in Breit Wigner block(s) " << endl; 
      exit(1);
    }

    fA = fLight.GetA();
    fmu = (uma*fA*fn)/(fA+fn);
    fR = fR0*(TMath::Power(fA,(1.0/3.0)) + TMath::Power(fn,(1.0/3.0)));
    frho_0 = Calculate_Rho_Var(fX0);
    fP0 = Penetration_Factor(frho_0);
    fSF0 = Shift_Factor(frho_0);

    fQ = fLight.Mass() + fn*fNeutron.Mass() - fHeavy.Mass();
    cout << "Q de reaction = " << fQ << endl;

    
  }//end reading token in input file
}//end ReadConfigurationFile 


double BreitWigner::Eval(const double & x_) const {
  
  if(x_<0) return 0.;
  
  double rho_x, W, dX;
  rho_x = Calculate_Rho_Var(x_);
  W = Width(rho_x);
  dX = Shift(rho_x);

  return fparam0*(W*fW0/4.0)/(TMath::Power(fX0+dX-x_,2.0) + TMath::Power(W/2.0, 2.0));
}

double BreitWigner::Eval_For_TF1(const double * x_, const double * /* params_ */) const
{
    return this->Eval(*x_);  
}


TF1 * BreitWigner::Make_TF1() const
{
  return new TF1("BW_TF1", this,&BreitWigner::Eval_For_TF1, 0.0, 10.0, 0);
}

TH1D * BreitWigner::Make_TH1() const
{
  TF1 *func =  new TF1("BW_TF1", this, &BreitWigner::Eval_For_TF1, 0.0, 10.0, 0);
  TH1D *Hist = dynamic_cast<TH1D*>(func->GetHistogram());
  
  TAxis *AxisX = Hist->GetXaxis();
  AxisX->Set( AxisX->GetNbins(), fQ + (AxisX->GetXmin()), fQ + (AxisX->GetXmax()) );
  
  return Hist;
}

double BreitWigner::Calculate_Rho_Var( double E) const {
  return sqrt(2.0*fmu*E)*fR/hc;
}


double BreitWigner::Width(double rho_x) const {
  return fW0*Width_Ratio(rho_x);
}

double BreitWigner::Shift(double rho_x) const {
  return fW0*Shift_Ratio(rho_x);
}


double BreitWigner::Width_Ratio(double rho_x) const {
  
  if(fell>=0 and !fWidthCst){
    double fP = Penetration_Factor(rho_x);
    return fP/fP0;
  }
  else return 1.0;
}

double BreitWigner::Shift_Ratio(double rho_x) const {
  if(fell>=0 and !fShiftZero){
    return (fSF0 - Shift_Factor(rho_x))/(2.0*fP0);
  }
  else return 0.0;
}

double BreitWigner::Shift_Factor(double rho) const {
  if(rho==0) rho = 1e-36; //to avoid singularities
  
  double J = std::cyl_bessel_j(fell+0.5, rho);
  double Jm = std::cyl_bessel_j(fell-0.5, rho);
  double Y  = std::cyl_neumann(fell+0.5, rho);
  double Ym = std::cyl_neumann(fell-0.5, rho);
  
  double sum2 = J*J + Y*Y ;
   
  return -(sum2*fell-(J*Jm+Y*Ym)*rho)/sum2;  
}


double BreitWigner::Penetration_Factor(double rho) const {
  if(rho==0) rho = 1e-36; // to avoid singularities
  
  Double_t J = std::cyl_bessel_j(fell+0.5,rho);
  Double_t Y = std::cyl_neumann(fell+0.5,rho);

  return  (2./TMath::Pi())/(TMath::Power(J,2)+TMath::Power(Y,2)); // in principle we must differiente rho / k but  R is the same  
}



void BreitWigner::Dump() const {
  //print
}



//--------------- /!\ For test only /!\----------------- //
/*
TF1 * BreitWigner::Make_TF1(const char * name_) const
{
  return new TF1(name_, this, &BreitWigner::Eval_For_TF1, 0.0, 10.0, 0);
}
*/


void BreitWigner::Save_TF1() const {

  TF1 *BWFunc =  new TF1("BWFunc", this, &BreitWigner::Eval_For_TF1, 0.0, 10.0, 0);
  
  TFile myfile("../../Projects/S034_audrey/output/BreitWigner/testBW.root","RECREATE");
 
  BWFunc->Write();
  myfile.Close();
  delete BWFunc;
}



void BreitWigner::Save_Histo() const {

  TF1 *fonc =  new TF1("fonc", this, &BreitWigner::Eval_For_TF1, 0.0, 10.0, 0);
  
  TFile myfile("../../Projects/S034_audrey/output/BreitWigner/testHist.root","RECREATE");
  
  TH1 *hist = fonc->GetHistogram();

  TAxis *a = hist->GetXaxis();
  a->Set( a->GetNbins(), fQ + (a->GetXmin()), fQ + (a->GetXmax()) );
  
  hist->Write();
  myfile.Close();
  delete fonc;
}

#endif
