#ifndef __SAMURAIIDEALDATA__
#define __SAMURAIIDEALDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author:    contact address:                                      *
 *                                                                           *
 * Creation Date  :                                                          *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include <vector>
#include <string>


#include "TObject.h"
using namespace std ;


class TSamuraiIdealData : public TObject {
 private:
    
   vector <short> Detector_Number;
   /* 
      -1 = Magnet
      0 = FDC0
      1 = FDC1
      2 = FDC2
   */
   vector <double> Dep_Energy;  //Energy deposited in the detector material
   vector <double> Brho;
   //Position
   vector <double> Pos_X; // Exit position X-axis
   vector <double> Pos_Y; // Exit position Y-axis
   vector <double> Pos_Z; // Exit position Z-axis
   //Momentum
   vector <double> Mom_Mag; //Exit momentum magnitude
   vector <double> Mom_Theta; //Exit momentum Theta
   vector <double> Mom_Phi; //Exit momentum Phi

   
 public:
   TSamuraiIdealData();
   virtual ~TSamuraiIdealData();

   void   Clear();
   void   Clear(const Option_t*) {};
   void   Dump() const;

   /////////////////////           GETTERS           ////////////////////////
   unsigned int   GetMult()            const {return Detector_Number.size();}
   short          GetDetNumber(int i)  const {return Detector_Number[i];}
   double         GetDepEnergy(int i)  const {return Dep_Energy[i];}
   double         GetBrho(int i)     const {return Brho[i];}
   //Position
   double         GetPosX(int i)       const {return Pos_X[i];}
   double         GetPosY(int i)       const {return Pos_Y[i];}
   double         GetPosZ(int i)       const {return Pos_Z[i];}
   //Momentum
   double         GetMomMag(int i)     const {return Mom_Mag[i];}
   double         GetMomTheta(int i)   const {return Mom_Theta[i];}
   double         GetMomPhi(int i)     const {return Mom_Phi[i];}


   /////////////////////           SETTERS           ////////////////////////
   //void SetData (short detector, double energy, G4ThreeVector pos, G4ThreeVector //mom, double brho){
   //   SetData(detector, energy, pos.x(), pos.y(), pos.z(), mom.getR(), mom.getTheta()//, mom.getPhi(), brho);
   //}
   void SetData (short detector, double energy, double pos_x, double pos_y, 
         double pos_z, double mom_r, double mom_theta, double mom_phi, double brho);

   void SetDetectorNumber(short i)  {Detector_Number.push_back(i);}
   void SetEnergy (double energy)   {Dep_Energy.push_back(energy);}
   ClassDef(TSamuraiIdealData,1)  // TSamuraiIdealData structure
};

#endif
