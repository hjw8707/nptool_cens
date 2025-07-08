/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: ValerianAlcindor  contact address:
 *valcindor@@ikp.tu-darmstadt.de
 *                                                                           *
 * Creation Date  : September 2021                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                  []                                             *
 *  A quite Standard Geant4 EventAction class.                               *
 *  Call the Fill method of the output tree.                                 *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "SteppingAction.hh"

#include "G4UnitsTable.hh"
#include "NPOptionManager.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction() { m_cut_parent_id = NPOptionManager::getInstance()->GetCutParentID(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    if (track->GetParentID() > m_cut_parent_id) track->SetTrackStatus(fKillTrackAndSecondaries);
}
