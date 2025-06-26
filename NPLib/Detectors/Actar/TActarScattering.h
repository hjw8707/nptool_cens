#ifndef TActarSCATTERING_H
#define TActarSCATTERING_H
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

class TActarScattering : public TObject {

public:
    TActarScattering();
    ~TActarScattering() {};

    //////////////////////////////////////////////////////////////
    // Inherited from TObject and overriden to avoid warnings
public:
    void Clear() {};
    void Clear(const Option_t*) {};

    //////////////////////////////////////////////////////////////
    // data obtained after reconstruction of an event

    double Theta;
    double Phi;
    double Range;
    double Charge;
    double Energy;
    TVector3 Vertex;
    TVector3 EndOfTrack;
    double ThetaShort;
    double PhiShort;

    int Voxels;
    int VoxelsInTrigger;

public:
    ClassDef(TActarScattering, 1) //ActarScattering structure
};


#endif