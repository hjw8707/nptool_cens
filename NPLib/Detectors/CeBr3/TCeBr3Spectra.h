#ifndef TCeBr3SPECTRA_H
#define TCeBr3SPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : December 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CeBr3 Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TCeBr3Data.h"
#include "TCeBr3Physics.h"

// Forward Declaration
class TCeBr3Physics;


class TCeBr3Spectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TCeBr3Spectra();
    TCeBr3Spectra(unsigned int NumberOfDetectors);
    ~TCeBr3Spectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TCeBr3Data*);
    void FillPreTreatedSpectra(TCeBr3Data*);
    void FillPhysicsSpectra(TCeBr3Physics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
