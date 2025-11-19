
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TTree.h"
void macro2() {
  TFile* file = new TFile("root/ana/sunji_aa.root");
  TTree* tree = (TTree*)file->Get("PhysicsTree");

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
  c1->Divide(2, 2);

  const Int_t nReaction = 2;

  Int_t nEachEvent = 10000;
  Int_t nStartEventNumber[nReaction] = {0, 10000};
  const char* reactionName[nReaction] = {"aa_gs", "aa_ex"};

  // Missing Mass
  c1->cd(1);
  for (Int_t i = 0; i < nReaction; i++) {
    tree->Draw(Form("RecoilExcitationEnergy>>hREE_%s(50,-5,5)", reactionName[i]), 0, "goff", nEachEvent,
               nStartEventNumber[i]);
    TH1F* hREE = (TH1F*)gDirectory->Get(Form("hREE_%s", reactionName[i]));
    hREE->SetLineColor(i + 1);
    hREE->SetMarkerColor(i + 1);
    hREE->SetMarkerStyle(20);
    hREE->SetMarkerSize(2);
    hREE->Draw("same");
  }

  // dE-E plot
  c1->cd(2);
  for (Int_t i = 0; i < nReaction; i++) {
    tree->Draw(Form("dE:E>>hDEE_%s(200,0,80,200,0,50)", reactionName[i]), 0, "goff", nEachEvent, nStartEventNumber[i]);
    TH2F* hDE = (TH2F*)gDirectory->Get(Form("hDEE_%s", reactionName[i]));
    hDE->SetMarkerColor(i + 1);
    hDE->SetMarkerStyle(20);
    hDE->SetMarkerSize(0.5);
    hDE->Draw("scat,same");
  }

  // theta vs E plot
  c1->cd(3);
  for (Int_t i = 0; i < nReaction; i++) {
    tree->Draw(Form("OutgoingEnergy:OutgoingThetaLab/3.1415926*180>>hThetaE_%s(360,0,180,160,0,80)", reactionName[i]),
               0, "goff", nEachEvent, nStartEventNumber[i]);
    TH2F* hThetaE = (TH2F*)gDirectory->Get(Form("hThetaE_%s", reactionName[i]));
    hThetaE->SetMarkerColor(i + 1);
    hThetaE->SetMarkerStyle(20);
    hThetaE->SetMarkerSize(0.1);
    hThetaE->Draw("scat,same");
  }

  // legend
  c1->cd(4);
  TLegend* legend = new TLegend(0.2, 0.2, 0.9, 0.9);
  for (Int_t i = 0; i < nReaction; i++) {
    TH1F* hREE = (TH1F*)gDirectory->Get(Form("hREE_%s", reactionName[i]));
    legend->AddEntry(hREE, reactionName[i], "lp");
  }
  legend->Draw();
}