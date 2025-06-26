#ifndef TACTARBEAM_H
#define TACTARBEAM_H
/*****************************************************************************
 * Copyright (C) 2009-2017   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Damien Thisse  contact address: damien.thisse@cea.fr     *
 *                                                                           *
 * Creation Date  : November 2024                                            *
 * Last update    : November 2024                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Actar Scattring events properties                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
// C++ headers

using namespace std;

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

class TActarBeam : public TObject {

public:
    TActarBeam();
    ~TActarBeam() {};

    //////////////////////////////////////////////////////////////
    // Inherited from TObject and overriden to avoid warnings
public:
    void Clear() {};
    void Clear(const Option_t*) {};

    void SetPoints(TVector3 & EntryPoint, TVector3 & ExitPoint); //!
    void SetVoxels(int Voxels); //!
    void SetIsScatteredBeam(bool IsScatteredBeam); //!
    //////////////////////////////////////////////////////////////
    // data obtained after reconstruction of an event

    TVector3 f_EntryPoint;
    TVector3 f_ExitPoint;

    double f_Theta;
    int f_Voxels;
    bool f_IsScatteredBeam;
    


public:
    ClassDef(TActarBeam, 1) //ActarBeam structure
};


#endif