#ifndef TCACAOSPECTRA_H
#define TCACAOSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: jwhwang@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : 4월 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CACAO Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TCACAOData.h"
#include "TCACAOPhysics.h"

// Forward Declaration
class TCACAOPhysics;


class TCACAOSpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TCACAOSpectra();
    TCACAOSpectra(unsigned int NumberOfDetectors);
    ~TCACAOSpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TCACAOData*);
    void FillPreTreatedSpectra(TCACAOData*);
    void FillPhysicsSpectra(TCACAOPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
