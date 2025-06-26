/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Damien THISSE contact address: damien.thisse@cea.fr    *
 *                                                                           *
 * Creation Date   : October 2024                                            *
 * Last update     : October 2024                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class deal with finding all the track event by event                *
 *****************************************************************************/

#include "NPRansacACTAR.h"
#include "NPLinearRansac3D.h"

using namespace NPL;
using namespace std;

/////////////////////////////////////////////////
RansacACTAR::RansacACTAR()
{
    fRANSACMaxIteration = 1000;
    fRANSACThreshold = 16;//100;
    fRANSACPointThreshold = 0.05;//0.07;
    fRANSACChargeThreshold = 0.05;//0.07;
    fRANSACDistance = 7;//7;
    fMAXBarycenterDistance = 10;
    fAngleMax = 2;//degree
    
    fNumberOfPadsX = 128;
    fNumberOfPadsY = 128;

    Rand = new TRandom3();
}

vector<LinearCluster3D> RansacACTAR::SimpleRansac(vector<double> & vX, vector<double> & vY, vector<double> & vZ, vector<double> & vQ){
    fTotalCharge = 0;
    fOriginalCloudSize = vX.size();
    vector<int> index;
    vector<LinearCluster3D> all_clusters;

    for(int QQ = 0; QQ < vQ.size(); QQ++){
        fTotalCharge+= vQ.at(QQ);
        index.push_back(QQ);
    }
    double RemainingCharge = fTotalCharge;
    int loop_counter = 0;
    // cout << "Total charge is " << fTotalCharge << " ---> Threshold to end RANSAC is thus " <<  fTotalCharge*fRANSACChargeThreshold << endl;
    // cout << "Total voxels is " << index.size() << " ---> Threshold to end RANSAC is thus " <<  fOriginalCloudSize*fRANSACPointThreshold << endl;
    do{
        loop_counter++;
        //cout << "Loop " << loop_counter << ": remaining charge is " << RemainingCharge << "; number of points remaining : " << index.size() << endl;
        LinearCluster3D foundCluster(&vX, &vY, &vZ, &vQ);
        int highest_inliers = 0;
        vector<int> index_to_remove;

        for(int i = 0; i < fRANSACMaxIteration; i++){
            LinearCluster3D aCluster(&vX, &vY, &vZ, &vQ);
            vector<int> pos_index;
            UInt_t rand1 = Rand->Integer(index.size());
            UInt_t rand2;
            do{
                rand2 = Rand->Integer(index.size());
            } while(rand2 == rand1);

            TVector3 pt1(vX[index[rand1]], vY[index[rand1]], vZ[index[rand1]]);
            TVector3 pt2(vX[index[rand2]], vY[index[rand2]], vZ[index[rand2]]);
            
            for(int j = 0; j < index.size(); j++){
                TVector3 pt(vX[index[j]], vY[index[j]], vZ[index[j]]);
                if(MinimumDistancePointLine(pt1, pt2, pt) < fRANSACDistance){
                    aCluster.AddIndex(index[j]);
                    pos_index.push_back(j);

                }
            }
            if(aCluster.size() > highest_inliers){
                foundCluster = aCluster;
                highest_inliers = aCluster.size();
                index_to_remove = pos_index;
            }
        }

        if(highest_inliers > fRANSACThreshold){
           all_clusters.push_back(foundCluster);
           for(int i = index_to_remove.size()-1; i >-1; i--){
                RemainingCharge-= vQ.at(index[index_to_remove[i]]);
                index.erase(index.begin()+index_to_remove.at(i));
           }
        //    RemainingCharge = 0;
        //    for(int QQ = 0; QQ < index.size(); QQ++){
        //         RemainingCharge+= vQ.at(index.at(QQ));
        //     }
        }
        else break;

    } while(RemainingCharge > fTotalCharge*fRANSACChargeThreshold && index.size() > fOriginalCloudSize*fRANSACPointThreshold);

    return all_clusters;
}

void RansacACTAR::ReadParameterValue(string filename)
{
    bool ReadingStatus = false;
    
    ifstream RansacParamFile;
    RansacParamFile.open(filename.c_str());
    
    if(!RansacParamFile.is_open()){
        cout << "No Paramter File found for Ransac parameters" << endl;
        cout << "Using the default parameter:" << endl;
        cout << "RANSAC MaxIteration= " << fRANSACMaxIteration << endl;
        cout << "RANSAC Threshold=" << fRANSACThreshold << endl;
        cout << "RANSAC ChargeThreshold= " << fRANSACChargeThreshold << endl;
        cout << "RANSAC Distance= " << fRANSACDistance << endl;
        cout << "RANSAC PointThreshold= " << fRANSACPointThreshold << endl;
        cout << "MAX Barycenter Distance= " << fMAXBarycenterDistance << endl;
        cout << "Angle Max to merge tracks= " << fAngleMax << endl;
    }
    else{
        string LineBuffer, whatToDo, DataBuffer;
        while(!RansacParamFile.eof()){
            getline(RansacParamFile,LineBuffer);
            string name = "ConfigRansac";
            if(LineBuffer.compare(0, name.length(), name) == 0){
                cout << endl;
                cout << "**** ConfigRansac found" << endl;
                ReadingStatus = true;
            }
            while(ReadingStatus){
                whatToDo="";
                RansacParamFile >> whatToDo;
                // Search for comment symbol (%)
                if (whatToDo.compare(0, 1, "%") == 0) {
                    RansacParamFile.ignore(numeric_limits<streamsize>::max(), '\n' );
                }
                else if (whatToDo.compare(0,19,"RANSACMaxIteration=") == 0) {
                    RansacParamFile >> DataBuffer;
                    fRANSACMaxIteration = atoi(DataBuffer.c_str());
                    cout << "/// RANSAC MaxIteration= " << " " << fRANSACMaxIteration << " ///" << endl;
                }
                else if (whatToDo.compare(0,15,"RANSACDistance=") == 0) {
                    RansacParamFile >> DataBuffer;
                    fRANSACDistance = atoi(DataBuffer.c_str());
                    cout << "/// RANSAC Distance= " << " " << fRANSACDistance << " ///" << endl;
                }
                else if (whatToDo.compare(0,22,"RANSACChargeThreshold=") == 0) {
                    RansacParamFile >> DataBuffer;
                    fRANSACChargeThreshold = atof(DataBuffer.c_str());
                    cout << "/// RANSAC ChargeThreshold= " << " " << fRANSACChargeThreshold << " ///" << endl;
                }
                else if (whatToDo.compare(0,21,"RANSACPointThreshold=") == 0) {
                    RansacParamFile >> DataBuffer;
                    fRANSACPointThreshold = atof(DataBuffer.c_str());
                    cout << "/// RANSAC PointThreshold= " << " " << fRANSACPointThreshold << " ///" << endl;
                }
                else if (whatToDo.compare(0,16,"RANSACThreshold=") == 0) {
                    RansacParamFile >> DataBuffer;
                    fRANSACThreshold = atoi(DataBuffer.c_str());
                    cout << "/// RANSAC Threshold= " << " " << fRANSACThreshold << " ///" << endl;
                }
                else if (whatToDo.compare(0,22,"MAXBarycenterDistance=") == 0) {
                    RansacParamFile >> DataBuffer;
                    fMAXBarycenterDistance = atoi(DataBuffer.c_str());
                    cout << "/// Max Barycenter Distance= " << " " << fMAXBarycenterDistance << " ///" << endl;
                }
                else if (whatToDo.compare(0,16,"AngleMaxToMerge=") == 0) {
                    RansacParamFile >> DataBuffer;
                    fAngleMax = atoi(DataBuffer.c_str());
                    cout << "/// Angle Max to merge tracks= " << " " << fAngleMax << " ///" << endl;
                }
                else{
                    ReadingStatus = false;
                }
            }
        }
    }
}