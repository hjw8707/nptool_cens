#ifndef SteppingAction_h
#define SteppingAction_h 1
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2021                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  A quite Standard Geant4 SteppingAction class.                            *
 *  Call the Fill method of the output tree.                                 *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// G4 header defining G4 types
#include "globals.hh"

// NPL
#include "TTrackInfo.h"

// G4 header
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4UserSteppingAction.hh"

// Root
#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SteppingAction : public G4UserSteppingAction {
   public:
    SteppingAction();
    ~SteppingAction() {};

   public:
    void UserSteppingAction(const G4Step* step);

   private:
    int m_cut_parent_id;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
