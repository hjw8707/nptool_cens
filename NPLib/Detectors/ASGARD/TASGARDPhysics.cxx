/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ASGARD treated data                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// STL
#include <stdlib.h>

#include <cmath>
#include <limits>
using namespace std;

#include "TASGARDPhysics.h"

//   NPL
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
using namespace NPUNITS;

//   ROOT
#include "TChain.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TASGARDPhysics)
    ///////////////////////////////////////////////////////////////////////////
    TASGARDPhysics::TASGARDPhysics() {
    m_EventData = new TASGARDData;
    m_PreTreatedData = new TASGARDData;
    m_EventPhysics = this;

    m_Beta = TVector3(0, 0, 0.125);  // fragment beta
};

/////////////////////////////////////////////////
void TASGARDPhysics::BuildPhysicalEvent() {
    PreTreat();
    // Addback Map
    unsigned int mysize = Gamma_Energy.size();
    for (unsigned int i = 0; i < 16; i++) {
        for (unsigned int g = 0; g < mysize; g++) {
            /////////////////////////////////////////////////////////////////
            // About core: there is no separate hit in the core
            // if (Clover_Number[g] == i + 1 && Segment_Number[g] == 0) {
            //     m_map_E[i] += Gamma_Energy[g];
            //     if (Gamma_Energy[g] > m_map_Core_MaxE[i]) {
            //         m_map_Core_MaxE[i] = Gamma_Energy[g];
            //         m_map_Core_Crystal[i] = Crystal_Number[g];
            //     }
            // }
            /////////////////////////////////////////////////////////////////
            if (Clover_Number[g] == i + 1 && Segment_Number[g] > 0 && Segment_Number[g] < 9) {
                m_map_E[i] += Gamma_Energy[g];
                if (Gamma_Energy[g] > m_map_Segment_MaxE[i]) {
                    m_map_Segment_MaxE[i] = Gamma_Energy[g];
                    m_map_Segment_Crystal[i] = Crystal_Number[g];
                    m_map_Segment[i] = Segment_Number[g];
                }
            }
        }
    }

    // Final Addback and Doppler Correction
    int zero = 0;
    for (int i = 0; i < 16; i++) {
        if (m_map_E.find(i) != m_map_E.end()) {
            int clover = i + 1;
            TVector3 Pos;
            if (m_map_Segment_MaxE[i] > 0)
                Pos = GetSegmentPosition(clover, m_map_Segment_Crystal[i], m_map_Segment[i]);
            else if (m_map_Core_MaxE[i] > 0)
                Pos = GetSegmentPosition(clover, m_map_Core_Crystal[i], zero);

            if (Pos.Mag() != 0) {
                double E = GetDopplerCorrectedEnergy(m_map_E[i], Pos, m_Beta);
                AddBack_DC.push_back(E);
                AddBack_E.push_back(m_map_E[i]);
                AddBack_Theta.push_back(Pos.Angle(m_Beta) * 180. / 3.141592653589793);
                AddBack_Clover.push_back(clover);
                if (m_map_Segment_MaxE[i] > 0) {
                    AddBack_Crystal.push_back(m_map_Segment_Crystal[i]);
                    AddBack_Segment.push_back(m_map_Segment[i]);
                } else {
                    AddBack_Crystal.push_back(m_map_Core_Crystal[i]);
                    AddBack_Segment.push_back(0);
                }
                AddBack_X.push_back(Pos.X());
                AddBack_Y.push_back(Pos.Y());
                AddBack_Z.push_back(Pos.Z());
            }
        }
    }
}

/////////////////////////////////////////////////
void TASGARDPhysics::PreTreat() {
    static CalibrationManager* cal = CalibrationManager::getInstance();
    static string name;
    unsigned int mysize = m_EventData->GetMultiplicityGe();
    double Eraw, Energy;
    int clover, crystal, segment;
    for (unsigned int i = 0; i < mysize; i++) {
        Eraw = m_EventData->GetGeEnergy(i);
        if (Eraw > 0) {
            clover = m_EventData->GetGeCloverNbr(i);
            crystal = m_EventData->GetGeCrystalNbr(i);
            segment = m_EventData->GetGeSegmentNbr(i);
            name = "ASGARD/D" + NPL::itoa(clover) + "_CRY" + NPL::itoa(crystal) + "_SEG" + NPL::itoa(segment) + "_E";
            Energy = cal->ApplyCalibration(name, Eraw);
            Gamma_Energy.push_back(Energy);
            Clover_Number.push_back(clover);
            Crystal_Number.push_back(crystal);
            Segment_Number.push_back(segment);
            Gamma_Time.push_back(m_EventData->GetGeTimeLED(i));

            // Doppler correction for each segment
            TVector3 v_pos(GetSegmentPosition(clover, crystal, segment));
            double E_dopp = GetDopplerCorrectedEnergy(Eraw, v_pos, m_Beta);
            Gamma_Dopp_Energy.push_back(E_dopp);

            // Look for Associate BGO
            // bool BGOcheck = false ;
            // for(unsigned j = 0 ;  j <  m_EventData->GetMultiplicityBGO() ; j++){

            //   if( m_EventData->GetBGOCloverNbr(j)== m_EventData->GetGeCloverNbr(i) && m_EventData->GetBGOEnergy(j)>20
            //   )
            //     BGOcheck = true ;
            // }
            // BGO.push_back(BGOcheck);
        }
    }
}

/////////////////////////////////////////////////
TVector3 TASGARDPhysics::GetPositionOfInteraction(unsigned int& i) {
    return GetSegmentPosition(Clover_Number[i], Crystal_Number[i], Segment_Number[i]);
}
/////////////////////////////////////////////////
// original energy, position, beta
double TASGARDPhysics::GetDopplerCorrectedEnergy(double& energy, TVector3 position, TVector3& beta) {
    // renorm pos vector
    position.SetMag(1);
    m_GammaLV.SetPx(energy * position.X());
    m_GammaLV.SetPy(energy * position.Y());
    m_GammaLV.SetPz(energy * position.Z());
    m_GammaLV.SetE(energy);
    m_GammaLV.Boost(-beta);
    return m_GammaLV.Energy();
}

/////////////////////////////////////////////////
// Add clover at the standard position of the array
// Take as argument the standard clover Id.
void TASGARDPhysics::AddCloverStandard(vector<int> CloverId) {
    for (unsigned int i = 0; i < CloverId.size(); i++) {
        if (CloverId[i] == 1) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(45 * deg);
            Pos.SetPhi(22.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 2) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(45 * deg);
            Pos.SetPhi(112.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 3) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(45 * deg);
            Pos.SetPhi(202.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 4) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(45 * deg);
            Pos.SetPhi(292.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 5) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(90 * deg);
            Pos.SetPhi(22.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 6) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(90 * deg);
            Pos.SetPhi(67.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 7) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(90 * deg);
            Pos.SetPhi(112.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 8) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(90 * deg);
            Pos.SetPhi(157.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 9) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(90 * deg);
            Pos.SetPhi(202.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 10) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(90 * deg);
            Pos.SetPhi(247.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 11) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(90 * deg);
            Pos.SetPhi(292.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 12) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(90 * deg);
            Pos.SetPhi(337.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 13) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(135 * deg);
            Pos.SetPhi(22.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 14) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(135 * deg);
            Pos.SetPhi(112.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 15) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(135 * deg);
            Pos.SetPhi(202.5 * deg);
            m_CloverPosition[ID] = Pos;
        }

        else if (CloverId[i] == 16) {
            TVector3 Pos(0, 0, 1);
            int ID = (CloverId[i]);
            Pos.SetMag(145 * mm);
            Pos.SetTheta(135 * deg);
            Pos.SetPhi(292.5 * deg);
            m_CloverPosition[ID] = Pos;
        }
    }

    return;
}
/////////////////////////////////////////////////
void TASGARDPhysics::AddClover(unsigned int ID, double R, double Theta, double Phi) {
    TVector3 Pos(0, 0, 1);
    Pos.SetTheta(Theta);
    Pos.SetPhi(Phi);
    Pos.SetMag(R);
    m_CloverPosition[ID] = Pos;
    return;
}
/////////////////////////////////////////////////
TVector3 TASGARDPhysics::GetCloverPosition(int& CloverNbr) { return m_CloverPosition[CloverNbr]; }
/////////////////////////////////////////////////
TVector3 TASGARDPhysics::GetCorePosition(int& CloverNbr, int& CoreNbr) {
    static double offset = 33.4;  // mm
    static double depth = 45;
    static TVector3 Pos;
    TVector3 CloverPos = m_CloverPosition[CloverNbr];

    if (CoreNbr == 0)
        Pos.SetXYZ(-offset, offset, depth);
    else if (CoreNbr == 1)
        Pos.SetXYZ(offset, offset, depth);
    else if (CoreNbr == 2)
        Pos.SetXYZ(offset, -offset, depth);
    else if (CoreNbr == 3)
        Pos.SetXYZ(-offset, -offset, depth);
    else
        cout << "Warning: ASGARD crystal number " << CoreNbr << " is out of range (1 to 4)" << endl;

    // Define reference axis as the clover direction
    Pos.RotateUz(CloverPos.Unit());
    Pos += CloverPos;
    return (Pos);
}
/////////////////////////////////////////////////
TVector3 TASGARDPhysics::GetSegmentPosition(int& CloverNbr, int& CoreNbr, int& SegmentNbr) {
    static double offsetXY1 = 10.4;  // mm
    static double offsetXY2 = 16.7;  // mm
    static double offsetZ1 = 15.5;   // mm
    static double offsetZ2 = 60.5;   // mm
    TVector3 CorePos = GetCorePosition(CloverNbr, CoreNbr);
    TVector3 CloverPos = GetCloverPosition(CloverNbr);
    static TVector3 Pos;

    if (SegmentNbr == 0 || SegmentNbr == 9)
        return GetCorePosition(CloverNbr, CoreNbr);
    else if (SegmentNbr == 1)
        Pos.SetXYZ(-offsetXY1, offsetXY1, offsetZ1);
    else if (SegmentNbr == 2)
        Pos.SetXYZ(offsetXY1, offsetXY1, offsetZ1);
    else if (SegmentNbr == 3)
        Pos.SetXYZ(offsetXY1, -offsetXY1, offsetZ1);
    else if (SegmentNbr == 4)
        Pos.SetXYZ(-offsetXY1, -offsetXY1, offsetZ1);
    else if (SegmentNbr == 5)
        Pos.SetXYZ(-offsetXY2, offsetXY2, offsetZ2);
    else if (SegmentNbr == 6)
        Pos.SetXYZ(offsetXY2, offsetXY2, offsetZ2);
    else if (SegmentNbr == 7)
        Pos.SetXYZ(offsetXY2, -offsetXY2, offsetZ2);
    else if (SegmentNbr == 8)
        Pos.SetXYZ(-offsetXY2, -offsetXY2, offsetZ2);
    else
        cout << "Warning: ASGARD segment number " << SegmentNbr << " is out of range (0 to 9)" << endl;

    // Each crystal is a rotation of the previous one
    if (CoreNbr == 2)
        Pos.RotateZ(90 * deg);
    else if (CoreNbr == 3)
        Pos.RotateZ(180 * deg);
    else if (CoreNbr == 4)
        Pos.RotateZ(270 * deg);

    // Define reference axis as the core position
    Pos.RotateUz(CorePos.Unit());
    Pos += CorePos;
    return (Pos);
}

/////////////////////////////////////////////////
void TASGARDPhysics::ReadConfiguration(NPL::InputParser parser) {
    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("ASGARD", "Clover");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " free clovers found " << endl;

    vector<string> token = {"CloverID", "R", "Theta", "Phi",
                            "Beta"};  // Beta not used,
                                      // but introduced for consistency with asgard-simu
    for (unsigned int i = 0; i < blocks.size(); i++) {
        if (blocks[i]->HasTokenList(token)) {
            double R = blocks[i]->GetDouble("R", "mm");
            double Theta = blocks[i]->GetDouble("Theta", "deg");
            double Phi = blocks[i]->GetDouble("Phi", "deg");
            vector<double> beta = blocks[i]->GetVectorDouble("Beta", "deg");
            int id = blocks[i]->GetInt("CloverID");
            AddClover(id, R, Theta, Phi);
        }

        else {
            cout << "Warning: check your input file formatting " << endl;
        }
    }

    blocks.clear();
    blocks = parser.GetAllBlocksWithTokenAndValue("ASGARD", "Standard");
    token.clear();
    token = {"CloverID"};

    for (unsigned int i = 0; i < blocks.size(); i++) {
        if (blocks[i]->HasTokenList(token)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << "//// Standard clovers found " << endl;
            vector<int> id = blocks[i]->GetVectorInt("CloverID");

            AddCloverStandard(id);
        }

        else {
            cout << "Warning: check your input file formatting " << endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////
void TASGARDPhysics::InitializeRootInputRaw() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("ASGARD", true);
    if (inputChain->FindBranch("fASG_*")) inputChain->SetBranchStatus("fASG_*", true);
    inputChain->SetBranchAddress("ASGARD", &m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TASGARDPhysics::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    outputTree->Branch("ASGARD", "TASGARDPhysics", &m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TASGARDPhysics::Clear() {
    Gamma_Energy.clear();
    Gamma_Dopp_Energy.clear();
    Gamma_Time.clear();
    Crystal_Number.clear();
    Clover_Number.clear();
    Segment_Number.clear();
    //  BGO.clear();

    AddBack_E.clear();
    AddBack_DC.clear();
    AddBack_Theta.clear();
    AddBack_Clover.clear();
    AddBack_Crystal.clear();
    AddBack_Segment.clear();
    AddBack_X.clear();
    AddBack_Y.clear();
    AddBack_Z.clear();

    m_map_E.clear();
    m_map_Core_Crystal.clear();
    m_map_Core_MaxE.clear();
    m_map_Segment_Crystal.clear();
    m_map_Segment.clear();
    m_map_Segment_MaxE.clear();
}
///////////////////////////////////////////////////////////////////////////
void TASGARDPhysics::ClearEventData() {
    m_EventData->Clear();
    m_PreTreatedData->Clear();
}
///////////////////////////////////////////////////////////////////////////
void TASGARDPhysics::AddParameterToCalibrationManager() {
    CalibrationManager* Cal = CalibrationManager::getInstance();
    for (int i = 0; i < 16; ++i) {
        for (int cry = 0; cry < 4; cry++) {
            // core are 0 and 9 , segment 1 to 8
            for (int j = 0; j < 10; ++j) {
                Cal->AddParameter(
                    "ASGARD", "D" + NPL::itoa(i + 1) + "_CRY" + NPL::itoa(cry + 1) + "_SEG" + NPL::itoa(j) + "_E",
                    "ASGARD_D" + NPL::itoa(i + 1) + "_CRY" + NPL::itoa(cry + 1) + "_SEG" + NPL::itoa(j) + "_E");
            }
        }
    }
    return;
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TASGARDPhysics::Construct() { return (NPL::VDetector*)new TASGARDPhysics(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_asgard {
   public:
    proxy_asgard() {
        NPL::DetectorFactory::getInstance()->AddToken("ASGARD", "ASGARD");
        NPL::DetectorFactory::getInstance()->AddDetector("ASGARD", TASGARDPhysics::Construct);
    }
};

proxy_asgard p;
}
