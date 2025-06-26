#ifndef PROFILE_EVAL
#define PROFILE_EVAL
/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: T.Efremov  contact address: theodore.efremov@cea.fr   *
 *                                                                           *
 * Creation Date  : 22/01/24                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Evaluator for 2D Profile                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include <TFile.h>
#include <TProfile2D.h>
#include <iostream>
using namespace std;
  
class ProfileEvaluator {
    private:
        TFile* fFile;
        TProfile2D* fProfile;

    public:

        //***********************************************************//

        ProfileEvaluator() : fFile(nullptr), fProfile(nullptr) {}

        //***********************************************************//

        bool LoadProfile(const char* filename, const char* profname) {
            fFile = TFile::Open(filename, "READ");
            if (!fFile || fFile->IsZombie()) return false;

            fProfile = (TProfile2D*)fFile->Get(profname)->Clone("Profile");
            fProfile->SetDirectory(0);
            return fProfile != nullptr;
        }

        //***********************************************************//

        Double_t Evaluate(Double_t x, Double_t y, Bool_t useInterpolation = true) {
            if (!fProfile) return 0.0;

            Int_t binx = fProfile->GetXaxis()->FindBin(x);
            Int_t biny = fProfile->GetYaxis()->FindBin(y);

            if (!useInterpolation) return fProfile->GetBinContent(binx, biny);

            //Def Bin dist
            Double_t BinDistX = fProfile->GetXaxis()->GetBinCenter(2) - fProfile->GetXaxis()->GetBinCenter(1);
            Double_t BinDistY = fProfile->GetYaxis()->GetBinCenter(2) - fProfile->GetYaxis()->GetBinCenter(1);

            // Normalisation var
            Double_t norma = 0;

            // Fetch Bin Center for surrounding 3 by 3 bin and value
            vector<Double_t> Cx(3) , Cy(3) ;
            vector<vector<Double_t>> Q(3, vector<Double_t>(3)), W(3,vector<Double_t>(3));        

            for (int i = -1 ; i<=1 ; i++){
                // Loading nearby bins
                int TempX = binx +i ;
                int PosX = i + 1 ;

                // Edge cases handling    
                if (TempX < 1) TempX = 1 ; 
                if (TempX > fProfile->GetNbinsX()) TempX = fProfile->GetNbinsX() ;

                // Loading center
                Cx.at(PosX) = fProfile->GetXaxis()->GetBinCenter(TempX);
                //Normat that we use POSX here

                // Loading value and weight
                for (int j = -1 ; j<=1 ; j++){
                    // Need to translate to acess the vector
                    int PosY = j + 1 ;
                    int TempY = biny + j ;

                    if (TempY < 1) TempY = 1 ;  
                    if (TempY > fProfile->GetNbinsY()) TempY = fProfile->GetNbinsY() ; 

                    Cy.at(PosY) = fProfile->GetYaxis()->GetBinCenter(TempY); 
                   
                    // Loading value
                    Q.at(PosX).at(PosY) = fProfile->GetBinContent(TempX,TempY);

                    // Calculating weight
                    Double_t Wx = 0 , Wy = 0 ;
                    Wx = 1 - ( abs( x - Cx.at(PosX)) / (BinDistX) );
                    Wy = 1 - ( abs( y - Cy.at(PosY)) / (BinDistY) );
                   
                    //handling max range 
                    if(Wx < 0 || Wy < 0) W.at(PosX).at(PosY) = 0;
                    else W.at(PosX).at(PosY) = Wx * Wy;
                    // Verifiy that sum of the weight is 1
                    norma += W.at(PosX).at(PosY) ;
                }
            } 

            // now just calculate the return value

            Double_t InterPol = 0 ;

            for (int PosX = 0 ; PosX <= 2 ; PosX ++ ){
                for (int PosY = 0 ; PosY <= 2 ; PosY ++ ){
                    InterPol += Q.at(PosX).at(PosY) * W.at(PosX).at(PosY);
                }
            }
            return InterPol;
        }

        //***********************************************************//

        void PrintInfo() {
            if (!fProfile) return;
            printf("Profile: %s\n", fProfile->GetName());
            printf("X bins: %d, range: [%.2f, %.2f]\n", 
                    fProfile->GetNbinsX(),
                    fProfile->GetXaxis()->GetXmin(),
                    fProfile->GetXaxis()->GetXmax());
            printf("Y bins: %d, range: [%.2f, %.2f]\n", 
                    fProfile->GetNbinsY(),
                    fProfile->GetYaxis()->GetXmin(),
                    fProfile->GetYaxis()->GetXmax());
        }

        //***********************************************************//

        TProfile2D* GetProfile() { return fProfile; }

        //***********************************************************//
        ~ProfileEvaluator() {
        }
};
#endif
