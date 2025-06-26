#include <iostream>
#include <ctime>
#include <cstdlib>
using namespace std;

// ROOT headers
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TEllipse.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

// nptool headers
#include "TReactionConditions.h"
#include "TInteractionCoordinates.h"
#include "NPReaction.h"
#include "NPQFS.h"

using namespace NPL;

TCanvas* canvas1;
TCanvas* canvas2;

void CheckSimu(const char * fname = "Example5"){
  // for the style 
  gStyle->SetOptStat(0);
  gROOT->SetStyle("nptool");     
  gROOT->ForceStyle(true);  
  gStyle->SetPalette(1);

  // Open output ROOT file from NPTool simulation run
  TString path = gSystem->Getenv("NPTOOL");
  path += "/Outputs/Simulation/";
  TString inFileName = fname;
  if (!inFileName.Contains("root")) inFileName += ".root";
  TFile *inFile = new TFile(path + inFileName);
  if (!inFile->IsOpen()) exit(1);
  TTree *tree   = (TTree*) inFile->Get("SimulatedTree");

  // Connect the branches of the TTree and activate then if necessary
  // TReactionConditions branch
  TReactionConditions *reacCond = 0;
  tree->SetBranchAddress("ReactionConditions", &reacCond);
  tree->SetBranchStatus("ReactionConditions", 1);
 
  // TInteractionCoordinates branch
  TInteractionCoordinates *interCoord = 0;
  tree->SetBranchAddress("InteractionCoordinates", &interCoord);
  tree->SetBranchStatus("InteractionCoordinates", 1);

  // Declare histograms
  // for emitted particle
  TH1F *hEmThetaCM = new TH1F("hEmThetaCM", "Light Ejectile Theta CM",180, 0, 180);
  TH1F *hEmPhiCM = new TH1F("hEmPhiCM", "Light Ejectile Phi CM",370, -185, 185);
  TH1F *hEmInternalMomX = new TH1F("hEmInternalMomX", "Internal Momentum (X) of removed cluster",500, -500, 500);
  TH1F *hEmInternalMomY = new TH1F("hEmInternalMomY", "Internal Momentum (Y) of removed cluster",500, -500, 500);
  TH1F *hEmInternalMomZ = new TH1F("hEmInternalMomZ", "Internal Momentum (Z) of removed cluster",500, -500, 500);

  TH1F *hEmTheta1 = new TH1F("hEmTheta1", "Ejectile1 Theta (reaction frame)", 180, 0, 180);
  TH1F *hEmTheta2 = new TH1F("hEmTheta2", "Ejectile2 Theta (reaction frame)", 180, 0, 180);
  TH2F *hEmE1Theta1 = new TH2F("hEmE1Theta1",  "Kinematics (1)", 150, 0, 90, 250, 0, 500);
  TH2F *hEmE2Theta2 = new TH2F("hEmE2Theta2",  "Kinematics (2)", 150, 0, 90, 250, 0, 500);
  TH2F *hEmE1VsE2 = new TH2F("hEmE1VsE2", " E1 VS E2 (reaction frame)", 300, 0, 300,300,0,300);
  TH2F *hEmTheta1VsTheta2 = new TH2F("hEmTheta1VsTheta2", "Theta1 VS Theta2 (reaction frame)", 360, 0, 90,360,0,90);
  TH2F *hEmPhi1VsPhi2 = new TH2F("hEmPhi1VsPhi2", "Phi1 VS Phi2 (reaction frame)", 360, -180, 180,360,-180,180);

  TH1F *hEmTheta1IF = new TH1F("hEmTheta1IF","Ejectile1 Theta (reac. frame)", 180, 0, 180);
  TH1F *hEmPhi1IF   = new TH1F("hEmPhi1IF","Ejectile1 Phi (reac. frame)",360, 0, 360);
  TH2F *hEmE1Theta1IF = new TH2F("hEmE1Theta1IF","Kinematics (1)",150, 0, 90, 250, 0, 500);
  TH1F *hEmTheta2IF = new TH1F("hEmTheta2IF","Ejectile2 Theta (reac. frame)", 180, 0, 180);
  TH1F *hEmPhi2IF   = new TH1F("hEmPhi2IF","Ejectile2 Phi (reac. frame)", 360, 0, 360);
  TH2F *hEmE2Theta2IF = new TH2F("hEmE1Theta1IF","Kinematics (2)",150, 0, 90, 250, 0, 500);

  TH1F *hEmTheta1WF = new TH1F("hEmTheta1WF","Ejectile1 Theta (world frame)", 180, 0, 180);
  TH1F *hEmPhi1WF   = new TH1F("hEmPhi1WF","Ejectile1 Phi (world frame)", 360, 0, 360);
  TH1F *hEmTheta2WF = new TH1F("hEmTheta2WF","Ejectile2 Theta (world frame)", 180, 0, 180);
  TH1F *hEmPhi2WF   = new TH1F("hEmPhi2WF",   "Ejectile2 Phi (world frame)", 360, 0, 360);
  TH2F *hEmTheta1WFVsTheta2WF = new TH2F("hEmTheta1WFVsTheta2WF", "Theta1WF VS Theta2WF)", 360, 0, 90,360,0,90);

  TH1F *hEmOpAngle = new TH1F("hEmOpAngle","Open. angle between (1) and (2)", 100, 0, 100);
  TH1F *hdPhi      = new TH1F("hdPhi","Phi angle diff. between (1) and (2)", 200, 0, 200);


  // Read the TTree
  Int_t nentries = tree->GetEntries();
  cout << endl << " TTree contains " << nentries << " events" << endl;

  for (Int_t i = 0; i < nentries; i++) {
 // for (Int_t i = 0; i < 10; i++) {
    if (i%10000 == 0 && i!=0)  {
      cout.precision(5);
      Double_t percent = (Double_t)i/nentries ;
      cout  << "\r Progression:  " << percent*100 << " %" << flush;
    }
    else if (i==nentries-1)  cout << "\r Progression:" << " 100%" << endl;

    // Get entry
    tree->GetEntry(i);

    // Fill histos
    // ejected particles
    hEmThetaCM  -> Fill(reacCond->GetThetaCM());
    //hEmPhiCM  -> Fill(reacCond->GetPhiCM());
    hEmInternalMomX  -> Fill(reacCond->GetInternalMomentum().X());
    hEmInternalMomY  -> Fill(reacCond->GetInternalMomentum().Y());
    hEmInternalMomZ  -> Fill(reacCond->GetInternalMomentum().Z());

    hEmTheta1   -> Fill(reacCond->GetTheta(0));
    hEmTheta2   -> Fill(reacCond->GetTheta(1));
    hEmTheta1VsTheta2   -> Fill(reacCond->GetTheta(1), reacCond->GetTheta(0));
    hEmE1Theta1  -> Fill(reacCond->GetTheta(0), reacCond->GetKineticEnergy(0));
    hEmE2Theta2  -> Fill(reacCond->GetTheta(1), reacCond->GetKineticEnergy(1));
    hEmPhi1VsPhi2-> Fill(reacCond->GetPhi(1), reacCond->GetPhi(0));
    hEmE1VsE2-> Fill(reacCond->GetKineticEnergy(1), reacCond->GetKineticEnergy(0));

    hEmTheta1IF -> Fill(reacCond->GetTheta_BeamFrame(0));
    hEmPhi1IF -> Fill(reacCond->GetPhi_BeamFrame(0));
    hEmE1Theta1IF  -> Fill(reacCond->GetTheta_BeamFrame(0), reacCond->GetKineticEnergy(0));
    hEmTheta2IF  -> Fill(reacCond->GetTheta_BeamFrame(1));
    hEmPhi2IF -> Fill(reacCond->GetPhi_BeamFrame(1));
    hEmE2Theta2IF  -> Fill(reacCond->GetTheta_BeamFrame(1), reacCond->GetKineticEnergy(1));

    hEmTheta1WF  -> Fill(reacCond->GetTheta_WorldFrame(0));
    hEmTheta2WF  -> Fill(reacCond->GetTheta_WorldFrame(1));
    hEmTheta1WFVsTheta2WF -> Fill(reacCond->GetTheta_WorldFrame(1), reacCond->GetTheta_WorldFrame(0));
    hEmPhi1WF  -> Fill(reacCond->GetPhi_WorldFrame(0));
    hEmPhi2WF  -> Fill(reacCond->GetPhi_WorldFrame(1));

/*
    double theta1 = reacCond->GetTheta_BeamFrame(0)*TMath::Pi()/180.;
    double theta2 = reacCond->GetTheta_BeamFrame(1)*TMath::Pi()/180.;
    double phi1 = reacCond->GetPhi_BeamFrame(0)*TMath::Pi()/180.;
    double phi2 = reacCond->GetPhiL_BeamFrame(1)*TMath::Pi()/180.;
    double Opang = acos( sin(theta1)*sin(theta2)*cos(phi1-phi2) +
                         cos(theta1)* cos(theta2) );                      
*/
    double theta1 = reacCond->GetTheta(0)*TMath::Pi()/180.;
    double theta2 = reacCond->GetTheta(1)*TMath::Pi()/180.;
    double phi1 = reacCond->GetPhi(0)*TMath::Pi()/180.;
    double phi2 = reacCond->GetPhi(1)*TMath::Pi()/180.;
    double Opang = acos( sin(theta1)*sin(theta2)*cos(phi1-phi2) +
                         cos(theta1)* cos(theta2) );                      
    hEmOpAngle->Fill(Opang*180./TMath::Pi());

    double df = fabs(phi1-phi2)*180./TMath::Pi();

    if(df>0 && df <= 180) hdPhi->Fill(df);
    else hdPhi->Fill(360 - df);

    //reacCond->Dump();

  }

  /////////////////////////////////////////////////////
  // All properties in reaction frame (beam along Z)///
  /////////////////////////////////////////////////////
  
  canvas1 = new TCanvas("canvas1", "Reaction frame (beam along Z)",1000,1000);
  canvas1->Divide(3,3);

  canvas1->cd(1);
  hEmThetaCM->SetXTitle("#theta_{c.m.}");
  hEmThetaCM->SetYTitle("counts / 1^{#circ}");
  hEmThetaCM->GetYaxis()->SetTitleOffset(1.18);
  //hEmThetaCM->GetYaxis()->SetRangeUser(0,500);
  hEmThetaCM->Draw();
  
  canvas1->cd(2);
  hEmE1Theta1IF->SetXTitle("#theta_{1}");
  hEmE1Theta1IF->SetYTitle("E_{1} (MeV)");
  hEmE1Theta1IF->Draw("colz");

  canvas1->cd(3);
  hdPhi->Draw();
  hdPhi->SetXTitle("#phi_{1} - #phi_{2} (deg)");
  hdPhi->SetYTitle("Counts");

  canvas1->cd(4);
  hEmTheta1VsTheta2->Draw("colz");
  hEmTheta1VsTheta2->SetXTitle("#theta_{1} (deg)");
  hEmTheta1VsTheta2->SetYTitle("#theta_{2} (deg)");
  NPL::QFS qfs;
  qfs.ReadConfigurationFile("Example5.reaction");
  qfs.SetMomentumSigma(0.);
  TGraph* Kine = qfs.GetTheta2VsTheta1(1);
  Kine->SetLineWidth(2);
  Kine->SetLineColor(kRed);
  Kine->SetLineStyle(2);
  Kine->Draw("csame");
  TGraph* Kine2 = qfs.GetTheta2VsTheta1(10);
  Kine2->SetLineColor(kRed);
  Kine2->SetMarkerColor(kRed);
  Kine2->SetMarkerStyle(20);
  Kine2->SetMarkerSize(1.3);
  Kine2->Draw("Psame");

  canvas1->cd(5);
  hEmPhi1VsPhi2->Draw("colz");
  hEmPhi1VsPhi2->SetXTitle("#phi_{1} (deg)");
  hEmPhi1VsPhi2->SetYTitle("#phi_{2} (deg)");
  TGraph* KinePhi = qfs.GetPhi2VsPhi1(1);
  KinePhi->SetMarkerColor(kRed);
  KinePhi->SetMarkerSize(0.4);
  KinePhi->Draw("Psame");

  canvas1->cd(6);
  hEmOpAngle->Draw();
  hEmOpAngle->SetXTitle("Opening angle (1-2)  (deg)");
  hEmOpAngle->SetYTitle("Counts");

  canvas1->cd(7);
  hEmInternalMomX->Draw();
  hEmInternalMomX->SetXTitle("P_{x} (MeV/c)");
  hEmInternalMomX->SetYTitle("Counts");

  canvas1->cd(8);
  hEmInternalMomY->Draw();
  hEmInternalMomY->SetXTitle("P_{y} (MeV/c)");
  hEmInternalMomY->SetYTitle("Counts");

  canvas1->cd(9);
  hEmInternalMomZ->Draw();
  hEmInternalMomZ->SetXTitle("P_{z} (MeV/c)");
  hEmInternalMomZ->SetYTitle("Counts");

  ////////////////////////////////////
  // Reaction frame VS World Frame ///
  ////////////////////////////////////

  canvas2 = new TCanvas("canvas2", "Comparison Reaction frame VS World frame",1000,1000);
  canvas2->Divide(2,2);

  canvas2->cd(1);
  hEmTheta1->Draw();
  hEmTheta1->SetLineColor(kBlue);
  hEmTheta1->SetFillColorAlpha(kBlue,0.35);
  hEmTheta1WF->Draw("same");
  hEmTheta1WF->SetLineColor(kRed);
  hEmTheta1WF->SetFillColorAlpha(kRed,0.35);
  hEmTheta1->SetXTitle("Theta (deg)");
  hEmTheta1->SetYTitle("Counts");
  hEmTheta1->GetYaxis()->SetTitleOffset(1.18);
  canvas2->cd(2);
  hEmTheta1WFVsTheta2WF->Draw("colz");
  hEmTheta1WFVsTheta2WF->SetXTitle("#theta_{1} (deg)");
  hEmTheta1WFVsTheta2WF->SetYTitle("#theta_{2} (deg)");
  Kine->Draw("csame");
  Kine2->Draw("Psame");




}


