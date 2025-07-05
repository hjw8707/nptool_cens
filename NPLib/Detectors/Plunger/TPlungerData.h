#ifndef __PlungerDATA__
#define __PlungerDATA__
/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com                        *
 *                                                                           *
 * Creation Date  : 7월 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Plunger Raw data                                    *
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
#include "TVector3.h"

class TPlungerData : public TObject {
    //////////////////////////////////////////////////////////////
    // data members are hold into vectors in order
    // to allow multiplicity treatment
   private:
    // Particle flying to the stopper
    vector<string> fPlunger_ParticleName;
    vector<double> fPlunger_Velocity;
    vector<double> fPlunger_KineticEnergy;
    vector<TVector3> fPlunger_Position;

    //////////////////////////////////////////////////////////////
    // Constructor and destructor
   public:
    TPlungerData();
    ~TPlungerData();

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
    // Particle flying to the stopper
    inline void SetPlunger(const string& ParticleName, const double& Velocity, const double& KineticEnergy,
                           const TVector3& Position) {
        fPlunger_ParticleName.push_back(ParticleName);
        fPlunger_Velocity.push_back(Velocity);
        fPlunger_KineticEnergy.push_back(KineticEnergy);
        fPlunger_Position.push_back(Position);
    };  //!

    //////////////////////    GETTERS    ////////////////////////
    // Particle flying to the stopper
    inline UShort_t GetMultPlunger() const { return fPlunger_ParticleName.size(); }
    inline string GetParticleName(const unsigned int& i) const { return fPlunger_ParticleName[i]; }    //!
    inline double GetVelocity(const unsigned int& i) const { return fPlunger_Velocity[i]; }            //!
    inline double GetKineticEnergy(const unsigned int& i) const { return fPlunger_KineticEnergy[i]; }  //!
    inline TVector3 GetPosition(const unsigned int& i) const { return fPlunger_Position[i]; }          //!

    //////////////////////////////////////////////////////////////
    // Required for ROOT dictionnary
    ClassDef(TPlungerData, 1)  // PlungerData structure
};

#endif
