/*****************************************************************************
 * Original Author: Jongwon Hwang    contact address: jwhwang@ibs.re.kr      *
 * Creation Date  : May 2023                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold STARK Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "TSTARKPhysics.h"

#include <stdlib.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include "NPDetectorFactory.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
#include "TChain.h"
#include "TRotation.h"

using namespace NPUNITS;
using namespace std;

ClassImp(TSTARKPhysics)

    ///////////////////////////////////////////////////////////////////////////
    TSTARKPhysics::TSTARKPhysics() {
    m_EventData = new TSTARKData;
    m_PreTreatedData = new TSTARKData;
    m_EventPhysics = this;

    ////////////////////////////////////////////////////////////
    // for geometry
    X6_activeX = 40.30;
    X6_activeY = 75.00;  // mm
    X6_NFrontStrips = 8;
    X6_NBackStrips = 4;

    BB10_activeX = 39.45;
    BB10_activeY = 74.15;  // mm
    BB10_NFrontStrips = 8;
    BB10_NBackStrips = 1;

    QQQ5_activeInR = 25.25;
    QQQ5_activeOutR = 81.95;  // mm
    QQQ5_NRStrips = 32;
    QQQ5_NAStrips = 4;
    ////////////////////////////////////////////////////////////
}

///////////////////////////////////////////////////////////////////////////
TSTARKPhysics::~TSTARKPhysics() {}

void TSTARKPhysics::AddDetector(string Type, TVector3 Pos, int Flip, int Rev, double Beta, int Group) {
    m_Type.push_back(Type);
    m_DetPos.push_back(Pos);
    m_Flip.push_back(Flip);
    m_Rev.push_back(Rev);
    m_Beta.push_back(Beta);
    m_Group.push_back(Group);
    m_GroupLayer.push_back(0);  // initilization with 0 (not in any group)
    ////////////////////////////////////////////////////////////
    // Detector orientation (only useful for X6, but......)
    TVector3 V(0, 1, 0);  // Y
    if (Type == "QQQ5") {
        V.RotateZ(Beta);
        if (Flip) V.RotateY(TMath::Pi());
    } else {  // for X6, BB10
        if (Flip) V.RotateY(TMath::Pi());
        if (Rev) V.RotateZ(TMath::Pi());
        V.RotateX(TMath::PiOver2());
        V.RotateZ(TMath::PiOver2() + Pos.Phi());
    }
    m_DetOri.push_back(V);
    ////////////////////////////////////////////////////////////

    AddStripPosition(Type, Pos, Flip, Rev, Beta);
}

void TSTARKPhysics::AddStripPosition(string Type, TVector3 Pos, int Flip, int Rev, double Beta) {
    //////////////////////////////////////////////////////
    // local coordinate, unit vector (rotation)
    TRotation rot;
    if (Type == "QQQ5") {
        rot.RotateZ(Beta);
        if (Flip) rot.RotateY(TMath::Pi());
    } else {  // for X6, BB10
        if (Flip) rot.RotateY(TMath::Pi());
        if (Rev) rot.RotateZ(TMath::Pi());
        rot.RotateX(TMath::PiOver2());
        rot.RotateZ(TMath::PiOver2() + Pos.Phi());
    }
    //////////////////////////////////////////////////////

    if (Type == "X6") {
        ////////////////////////////////////////////////////////////
        // Strip position calculation for X6
        ////////////////////////////////////////////////////////////
        double X6_frontStripPitch = X6_activeX / X6_NFrontStrips;
        double X6_backStripPitch = X6_activeY / X6_NBackStrips;

        vector<vector<TVector3>> frontStripPos;
        for (int i = 0; i < X6_NFrontStrips; i++) {
            vector<TVector3> backStripPos;
            for (int j = 0; j < X6_NBackStrips; j++) {
                TVector3 PosStrip;
                PosStrip.SetX(-X6_frontStripPitch * (i + 0.5) + X6_activeX / 2);
                PosStrip.SetY(X6_backStripPitch * (j + 0.5) - X6_activeY / 2);
                PosStrip = rot * PosStrip;
                PosStrip += Pos;
                backStripPos.push_back(PosStrip);
            }
            frontStripPos.push_back(backStripPos);
        }
        m_StripPos.push_back(frontStripPos);
        ////////////////////////////////////////////////////////////

    } else if (Type == "BB10") {
        ////////////////////////////////////////////////////////////
        // Strip position calculation for BB10
        ////////////////////////////////////////////////////////////
        double BB10_frontStripPitch = BB10_activeX / BB10_NFrontStrips;
        double BB10_backStripPitch = BB10_activeY / BB10_NBackStrips;

        vector<vector<TVector3>> frontStripPos;
        for (int i = 0; i < BB10_NFrontStrips; i++) {
            vector<TVector3> backStripPos;
            for (int j = 0; j < BB10_NBackStrips; j++) {
                TVector3 PosStrip;
                PosStrip.SetX(-BB10_frontStripPitch * (i + 0.5) + BB10_activeX / 2);
                PosStrip.SetY(BB10_backStripPitch * (j + 0.5) - BB10_activeY / 2);
                PosStrip = rot * PosStrip;
                PosStrip += Pos;
                backStripPos.push_back(PosStrip);
            }
            frontStripPos.push_back(backStripPos);
        }
        m_StripPos.push_back(frontStripPos);
        ////////////////////////////////////////////////////////////
    } else if (Type == "QQQ5") {
        ////////////////////////////////////////////////////////////
        // Strip position calculation for QQQ5
        ////////////////////////////////////////////////////////////
        double QQQ5_rStripPitch = (QQQ5_activeOutR - QQQ5_activeInR) / QQQ5_NRStrips;
        double QQQ5_aStripPitch = TMath::PiOver2() / QQQ5_NAStrips;

        vector<vector<TVector3>> frontStripPos;  // Radial Strip
        for (int i = 0; i < QQQ5_NRStrips; i++) {
            vector<TVector3> backStripPos;  // Annular strips
            for (int j = 0; j < QQQ5_NAStrips; j++) {
                Double_t r = QQQ5_activeOutR - QQQ5_rStripPitch * (i + 0.5);
                Double_t a = QQQ5_aStripPitch * (i + 0.5);
                TVector3 PosStrip;
                PosStrip.SetMagThetaPhi(r, TMath::PiOver2(), a);
                PosStrip = rot * PosStrip;
                PosStrip += Pos;
                backStripPos.push_back(PosStrip);
            }
            frontStripPos.push_back(backStripPos);
        }
        m_StripPos.push_back(frontStripPos);
        ////////////////////////////////////////////////////////////
    }

    ////////////////////////////////////////////////////////////
    // for checking
    //  auto lastDet = m_StripPos.back();
    //  for (auto it2 = lastDet.begin() ;
    //       it2 != lastDet.end() ; it2++) {
    //    size_t frontIdx = it2 - lastDet.begin();
    //    for (auto it3 = it2->begin() ;
    //	 it3 != it2->end() ; it3++) {
    //      size_t backIdx = it3 - it2->begin();
    //      std::cout << "(" << frontIdx << ", " << backIdx << "): ";
    //      it3->Print();
    //    }}
    ////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// BuildGroup()
//
// Define the layer order in the group (define m_GroupLayer)
// This function should be called after all detectors are added at the end of the ReadConfiguration()
///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::BuildGroup() {
    ////////////////////////////////////////////////////////////
    // m_GroupLayer[i] is already initialized with 0 (not in any group)
    ////////////////////////////////////////////////////////////

    // 1. find the number of groups
    int nGroup = *max_element(m_Group.begin(), m_Group.end());

    for (int i = 0; i < nGroup; i++) {
        vector<pair<double, int>> dist_det_idx;
        int iGroup = i + 1;  // group number from 1

        // 2. find the detectors belong to the group
        for (int j = 0; j < m_Group.size(); j++) {
            if (m_Group[j] == iGroup) {
                dist_det_idx.push_back(make_pair(m_DetPos[j].Mag(), j));
            }
        }

        // 3. find the layer order in the group (the layer near to the origin is the first)
        sort(dist_det_idx.begin(), dist_det_idx.end());
        for (size_t order = 0; order < dist_det_idx.size(); ++order) {
            int detIdx = dist_det_idx[order].second;
            m_GroupLayer[detIdx] = order + 1;  // order는 1부터 시작
        }
    }

    ////////////////////////////////////////////////////////////
    // for checking
    // for (int i = 0; i < m_Group.size(); i++) {
    //     cout << "detector " << i << " is in group " << m_Group[i] << " and layer " << m_GroupLayer[i] << endl;
    // }
    // ////////////////////////////////////////////////////////////
}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::BuildPhysicalEvent() {
    Clear();

    ////////////////////////////////////////////////////////////
    // for dE-E analysis
    std::map<int, std::map<int, double>> dE_map;  // (group number, (layer number, dE))
    std::map<int, double> E_map;                  // (group number, E)
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // Loop for All
    nhit = m_EventData->GetMult();
    for (Int_t i = 0; i < m_EventData->GetMult(); i++) {
        type[i] = m_EventData->GetType(i);  // 0: X6, 1: BB10, 2: QQQ5
        detN[i] = m_EventData->GetDetN(i);
        fStrN[i] = m_EventData->GetFStN(i);
        bStrN[i] = m_EventData->GetBStN(i);  // only 1 for BB10
        uppE[i] = m_EventData->GetUpE(i);
        dwnE[i] = m_EventData->GetDwE(i);
        sumE[i] = m_EventData->GetFrE(i);

        sPosArr.push_back(
            m_StripPos[detN[i] - 1][fStrN[i] - 1][bStrN[i] - 1]);  // detector and strip number from 1 (not 0)

        if (type[i] == 0) {  // for X6
            ////////////////////////////////////////
            // hit position with resistive strip
            //
            // 1. get the distance from the resistive strip center
            Double_t distOnStrip = (uppE[i] - dwnE[i]) / sumE[i] * X6_activeY / 2.;
            // 2. calculate strip center position
            // = average(strip center of ohmic #2 & strip center of ohmic #3)
            TVector3 stripCenter = m_StripPos[detN[i] - 1][fStrN[i] - 1][1];
            stripCenter += m_StripPos[detN[i] - 1][fStrN[i] - 1][2];
            stripCenter *= 0.5;
            // 3. final position = stripCenter + detOrientation(Y-axis)*distOnStrip
            TVector3 hPos = m_DetOri[detN[i] - 1];
            hPos *= distOnStrip;
            hPos += stripCenter;
            hPosArr.push_back(hPos);
            ////////////////////////////////////////
        } else {  // for BB10 & QQQ5
            hPosArr.push_back(
                m_StripPos[detN[i] - 1][fStrN[i] - 1][bStrN[i] - 1]);  // detector and strip number from 1 (not 0)
        }

        ////////////////////////////////////////////////////////////
        // for dE-E analysis
        int group = m_Group[detN[i] - 1];
        if (group == 0) continue;

        dE_map[group][m_GroupLayer[detN[i] - 1]] = sumE[i];
        E_map[group] += sumE[i];
        ////////////////////////////////////////////////////////////
    }
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // for dE-E analysis
    //

    // sort the group number by the total energy
    vector<pair<int, double>> group_E;
    for (auto& group : E_map) {
        group_E.push_back(make_pair(group.first, group.second));
    }
    sort(group_E.begin(), group_E.end(), [](const auto& a, const auto& b) { return a.second > b.second; });

    nGroup = group_E.size();
    for (int i = 0; i < nGroup; i++) {
        groupN[i] = group_E[i].first;
        groupE[i] = group_E[i].second;
    }

    for (int i = 0; i < min(20, nGroup); i++) {
        for (int j = 0; j < 20; j++) {
            groupdE[i][j] = dE_map[groupN[i]][j];
        }
    }
    ////////////////////////////////////////////////////////////

    return;
}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::Clear() {
    nhit = 0;
    sPosArr.clear();
    hPosArr.clear();
    nGroup = 0;
    for (int i = 0; i < 20; i++) {
        groupN[i] = 0;
        groupE[i] = 0;
        for (int j = 0; j < 20; j++) {
            groupdE[i][j] = 0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::ReadConfiguration(NPL::InputParser parser) {
    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("STARK");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    vector<string> reso = {"Type", "Reso"};
    vector<string> cart = {"Type", "POS"};
    vector<string> sphe = {"Type", "R", "Theta", "Phi"};
    vector<string> cyld = {"Type", "Rho", "Phi", "Z"};

    for (unsigned int i = 0; i < blocks.size(); i++) {
        ////////////////////////////////////////////////////////////
        // skip the resolution block
        if (blocks[i]->HasTokenList(reso)) continue;
        ////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////
        // Cartesian coordinate
        if (blocks[i]->HasTokenList(cart)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  STARK " << i + 1 << endl;
            string Type = blocks[i]->GetString("Type");
            TVector3 Pos = blocks[i]->GetTVector3("POS", "mm");
            int Flip = 0;
            if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
            int Rev = 0;
            if (blocks[i]->HasToken("Rev")) Rev = blocks[i]->GetInt("Rev");
            double Beta = 0;
            if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");
            int Group = 0;
            if (blocks[i]->HasToken("Group")) Group = blocks[i]->GetInt("Group");
            AddDetector(Type, Pos, Flip, Rev, Beta, Group);
        }
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Spherical coordinate
        else if (blocks[i]->HasTokenList(sphe)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  STARK " << i + 1 << endl;
            string Type = blocks[i]->GetString("Type");
            double R = blocks[i]->GetDouble("R", "mm");
            double Theta = blocks[i]->GetDouble("Theta", "deg");
            double Phi = blocks[i]->GetDouble("Phi", "deg");
            int Flip = 0;
            if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
            int Rev = 0;
            if (blocks[i]->HasToken("Rev")) Rev = blocks[i]->GetInt("Rev");
            double Beta = 0;
            if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");
            int Group = 0;
            if (blocks[i]->HasToken("Group")) Group = blocks[i]->GetInt("Group");
            TVector3 Pos;
            Pos.SetMagThetaPhi(R, Theta, Phi);
            AddDetector(Type, Pos, Flip, Rev, Beta, Group);
        }
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Cylindrical coordinate
        else if (blocks[i]->HasTokenList(cyld)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  STARK " << i + 1 << endl;
            string Type = blocks[i]->GetString("Type");
            double Rho = blocks[i]->GetDouble("Rho", "mm");
            double Phi = blocks[i]->GetDouble("Phi", "deg");
            double Z = blocks[i]->GetDouble("Z", "mm");
            int Flip = 0;
            if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
            int Rev = 0;
            if (blocks[i]->HasToken("Rev")) Rev = blocks[i]->GetInt("Rev");
            double Beta = 0;
            if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");
            int Group = 0;
            if (blocks[i]->HasToken("Group")) Group = blocks[i]->GetInt("Group");
            TVector3 Pos;
            Pos.SetXYZ(Rho * cos(Phi), Rho * sin(Phi), Z);
            AddDetector(Type, Pos, Flip, Rev, Beta, Group);
        }
        ////////////////////////////////////////////////////////////

        else {
            cout << "Error: check your input file formatting" << endl;
            exit(1);
        }
    }
    std::cout << "read complete" << std::endl;

    /////////////////////////////////////////////////////////////////////////////////
    // BuildGroup() should be called after all detectors are added
    BuildGroup();
    /////////////////////////////////////////////////////////////////////////////////
}

void TSTARKPhysics::AddParameterToCalibrationManager() {}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::InitSpectra() {}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::FillSpectra() {}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::CheckSpectra() {}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::ClearSpectra() {
    // To be done
}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::WriteSpectra() {}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::InitializeRootInputRaw() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("STARK", true);
    inputChain->SetBranchAddress("STARK", &m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::InitializeRootInputPhysics() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("STARK", true);
    inputChain->SetBranchAddress("STARK", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TSTARKPhysics::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    TFile* outputFile = RootOutput::getInstance()->GetFile();
    //    outputTree->Branch("STARK", "TSTARKPhysics", &m_EventPhysics);
    outputFile->WriteObjectAny(m_EventPhysics, "TSTARKPhysics", "STARK");

    outputTree->Branch("nhit", &nhit, "nhit/I");
    outputTree->Branch("type", type, "type[nhit]/I");
    outputTree->Branch("detN", detN, "detN[nhit]/I");
    outputTree->Branch("fStrN", fStrN, "fStrN[nhit]/I");
    outputTree->Branch("bStrN", bStrN, "bStrN[nhit]/I");
    outputTree->Branch("uppE", uppE, "uppE[nhit]/D");
    outputTree->Branch("dwnE", dwnE, "dwnE[nhit]/D");
    outputTree->Branch("sumE", sumE, "sumE[nhit]/D");
    outputTree->Branch("sPos", &sPosArr);
    outputTree->Branch("hPos", &hPosArr);

    outputTree->Branch("nGroup", &nGroup, "nGroup/I");
    outputTree->Branch("groupN", groupN, "groupN[nGroup]/I");
    outputTree->Branch("groupE", groupE, "groupE[nGroup]/D");
    outputTree->Branch("groupdE", groupdE, "groupdE[20][20]/D");  // should use fixed-length array for 2-dim array
}

////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSTARKPhysics::Construct() {
    // construct method to be passed to the DetectorFactory
    return (NPL::VDetector*)new TSTARKPhysics();
}

////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_STARK {
    // register the construct method to the factory
   public:
    proxy_STARK() {
        NPL::DetectorFactory::getInstance()->AddToken("STARK", "STARK");
        NPL::DetectorFactory::getInstance()->AddDetector("STARK", TSTARKPhysics::Construct);
    }
};

proxy_STARK p_STARK;
}
