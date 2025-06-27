#ifndef ISSArray_h
#define ISSArray_h 1
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Author: M. Labiche                     address: marc.labiche@stfc.ac.uk   *
 *                                                                           *
 * Creation Date  : July 2019                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This file describe the Iss charge particle Detector                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * Iss is a modular array made of DSSSD (Telescope). Each                    *
 *  Telescope is made of a single Stage:                                     *
 *  - A 300um Silicium, double-sided strip                                   *
 *  - possibility to add a second layer in future			     *
 *****************************************************************************/
#include "NPSVDetector.hh"
#include "TIssData.h"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "NPInputParser.h"
#include "G4GDMLParser.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ISS
{
   // Resolution
   const G4double ResoTime = 0.212765957    ;// = 500ps same as MUST for now      //   Unit is  ns/2.35
   //const G4double ResoTimeMust = 0.212765957    ;// = 500ps                 //   Unit is  ns/2.35
   const G4double ResoStrip    = 0.0149         ;// 0.0223 = 52keV of Resolution   //   Unit is MeV/2.35  14.861996
   const G4double TimeOffset   = 500            ;// 500 ns stop
   // Threshold
   const G4double ThresholdSi   = 50 * keV;
    // Geometry mother volume
   const G4double FaceWidth = 23.*mm ;
   const G4double FaceLength = 126*mm ;
   const G4double Length  = 1.008*mm ;

   const G4double AluStripThickness = 0.4*micrometer ;
   const G4double SiliconThickness  = 1000.*micrometer ;
   const G4double SiliconFaceWidth       = 22.*mm ;
   const G4double SiliconFaceLength      = 125.*mm ;
  const G4double VacBoxThickness   = 1.008*mm ;

   const G4int    Silicon_Front_NumberOfStrip = 128 ;
   const G4int    Silicon_Back_NumberOfStrip = 11 ;
   
   // Starting at the front and going to CsI
   const G4double AluStripFront_PosZ   = Length* -0.5 + 0.5*AluStripThickness;
   const G4double Silicon_PosZ      = AluStripFront_PosZ + 0.5*AluStripThickness + 0.5*SiliconThickness;
   const G4double AluStripBack_PosZ = Silicon_PosZ + 0.5*SiliconThickness + 0.5*AluStripThickness;
   const G4double VacBox_PosZ    = AluStripBack_PosZ + 0.5*AluStripThickness + 0.5* VacBoxThickness;
}

class IssArray : public NPS::VDetector
{
   ////////////////////////////////////////////////////
   /////// Default Constructor and Destructor /////////
   ////////////////////////////////////////////////////
public:
   IssArray()   ;
   virtual ~IssArray()   ;

   ////////////////////////////////////////////////////
   //////// Specific Function of this Class ///////////
   ////////////////////////////////////////////////////
public:
   // By Position Method
   void AddTelescope(   G4ThreeVector  TL       ,
                        G4ThreeVector  BL       ,
                        G4ThreeVector  BR       ,
                        G4ThreeVector  CT       ,
                        bool           wSi      );
   // By Angle Method
   void AddTelescope(   G4double    R        ,
                        G4double    Theta    ,
                        G4double    Phi      ,
                        G4double    beta_u   ,
                        G4double    beta_v   ,
                        G4double    beta_w   ,
                        bool        wSi      );

   // Effectively construct Volume
   // Avoid to have two time same code for Angle and Point definition
   void VolumeMaker( G4int             TelescopeNumber   ,
		     G4double 	       DetTgt		 ,
                     G4ThreeVector     MMpos             ,
                     G4RotationMatrix* MMrot             ,
                     bool              wSi               ,
                     G4LogicalVolume*  world             );


   ////////////////////////////////////////////////////
   /////////  Inherite from NPS::VDetector class ///////////
   ////////////////////////////////////////////////////
public:
   // Read stream at Configfile to pick-up parameters of detector (Position,...)
   // Called in DetecorConstruction::ReadDetextorConfiguration Method
   void ReadConfiguration(NPL::InputParser);

   // Construct detector and inialise sensitive part.
   // Called After DetecorConstruction::AddDetector Method
   void ConstructDetector(G4LogicalVolume* world);

   // Add Detector branch to the EventTree.
   // Called After DetecorConstruction::AddDetector Method
   void InitializeRootOutput();

   // Read sensitive part and fill the Root tree.
   // Called at in the EventAction::EndOfEventAvtion
   void ReadSensitive(const G4Event* event);


   ////////////////////////////////////////////////////
   ///////////Event class to store Data////////////////
   ////////////////////////////////////////////////////
private:
   TIssData* m_Event;

   ////////////////////////////////////////////////////
   ///////////////Private intern Data//////////////////
   ////////////////////////////////////////////////////
private:

  // For gdml geometries
   G4GDMLParser m_gdmlparser;

   // True if Define by Position, False is Define by angle
   vector<bool>   m_DefinitionType  ;

   // Used for "By Point Definition"
   vector<G4ThreeVector>   m_X1_Y1     ; // Top Left Corner Position Vector
   vector<G4ThreeVector>   m_X1_Y128   ; // Bottom Left Corner Position Vector
   vector<G4ThreeVector>   m_X128_Y1   ; // Bottom Right Corner Position Vector
   vector<G4ThreeVector>   m_X128_Y128 ; // Center Corner Position Vector

   // Used for "By Angle Definition"
   vector<G4double>  m_R      ; //  |
   vector<G4double>  m_Theta  ; //  > Spherical coordinate of Strips Silicium Plate
   vector<G4double>  m_Phi    ; //  |

   vector<G4double>  m_beta_u ; //  |
   vector<G4double>  m_beta_v ; //  > Tilt angle of the Telescope
   vector<G4double>  m_beta_w ; //  |

   // If Set to true if you want this stage on you telescope
   vector<bool>      m_wSi    ; // Silicium Strip 300um 128*128 Strip

   vector<bool>      m_wAddSi ; // Additionnal Thin Silicium Strip

   // Set to true if you want to see Telescope Frame in your visualisation
   bool m_non_sensitive_part_visiualisation ;
   
   
   ////////////////////////////////////////////////////
   ///////////////////// Scorer ///////////////////////
   ////////////////////////////////////////////////////
private:
   //   Initialize all Scorer used by the IssArray
   void InitializeScorers() ;

   //   Silicon Associate Scorer
   G4MultiFunctionalDetector* m_BOXScorer ;
   //G4MultiFunctionalDetector* m_StripScorer ;
   
    
    
   ////////////////////////////////////////////////////
   //////////////////// Material //////////////////////
   ////////////////////////////////////////////////////
private:
   //   Declare all material used by the IssArray
   void InitializeMaterial() ;
   // Si
   G4Material* m_MaterialSilicon;
   // Al
   G4Material* m_MaterialAluminium;
   // Iron
   G4Material* m_MaterialIron;
   // CsI
   G4Material* m_MaterialCsI;
   //  Vacuum
   G4Material* m_MaterialVacuum ;
   //  Mylar
   G4Material* m_MaterialMyl;
public:
    static NPS::VDetector* Construct();
};

extern G4RotationMatrix* Rotation(double tetaX, double tetaY, double tetaZ);
#endif
