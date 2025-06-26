#ifndef TSEASONSPECTRA_H
#define TSEASONSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Emmanuel Rey-herme                                       *
 * contact address: marine.vandebrouck@cea.fr                                *
 *                                                                           *
 * Creation Date  : septembre 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SEASON Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
 
 // NPLib headers
 #include "NPVSpectra.h"
 #include "TSEASONData.h"
 #include "TSEASONPhysics.h"

 // Forward Declaration
 class TSEASONPhysics;


 class TSEASONSpectra : public VSpectra {
   //////////////////////////////////////////////////////////////
   // constructor and destructor
   public:
     TSEASONSpectra();
     TSEASONSpectra(unsigned int NumberOfDetectors, unsigned int NumberOfStripsX, unsigned int NumberOfStripsY);
     ~TSEASONSpectra();

   //////////////////////////////////////////////////////////////
   // Initialization methods
   private:
     void InitRawSpectra();
     void InitPreTreatedSpectra();
     void InitPhysicsSpectra();

   //////////////////////////////////////////////////////////////
   // Filling methods
   public:
     void FillRawSpectra(TSEASONData*);
     void FillPreTreatedSpectra(TSEASONData*);
     void FillPhysicsSpectra(TSEASONPhysics*);

   //////////////////////////////////////////////////////////////
   // Detector parameters 
   private:
     unsigned int fNumberOfDetectors;
     unsigned int fNumberOfStripsX;
     unsigned int fNumberOfStripsY;
 };

 #endif
