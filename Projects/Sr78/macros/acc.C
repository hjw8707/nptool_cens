#include <iostream>
using namespace std;
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TText.h"
#include "TTree.h"
#include "TVector3.h"

TFile *dataFile = NULL;
TTree *dataTree = NULL;

TH1 *accHist = NULL;
TH1 *accTTTHist = NULL;
TH1 *accYY1Hist = NULL;
TH1 *calcHist = NULL;

void Merge(const char *type, double factor = 1) {
  const char *nn[4] = {"1st0", "1st2", "2nd0", "2nd2"};
  Int_t entries[4] = {40000, 8200, 1100, 660};
  // Int_t entries[4] = {100000, 20000, 2800, 1600};
  TFile *file[4];
  TFile *fileTemp[4];
  TTree *tree[4];
  TTree *treeTemp[4];

  // Cloning trees
  for (Int_t i = 0; i < 4; i++) {
    file[i] = new TFile(Form("root/ana/Sr78_ana_%s_%s.root", nn[i], type));
    tree[i] = static_cast<TTree *>(file[i]->Get("PhysicsTree"));
    fileTemp[i] = new TFile(Form("root/ana/Sr78_ana_%s_%s_temp.root", nn[i], type), "RECREATE");
    treeTemp[i] = tree[i]->CloneTree(int(entries[i] * factor));
    fileTemp[i]->Write();
    file[i]->Close();
    fileTemp[i]->Close();
    delete file[i];
    delete fileTemp[i];
    std::cout << "tree copied" << std::endl;
  }

  TList *list = new TList;
  for (Int_t i = 0; i < 4; i++) {
    file[i] = new TFile(Form("root/ana/Sr78_ana_%s_%s_temp.root", nn[i], type));
    tree[i] = static_cast<TTree *>(file[i]->Get("PhysicsTree"));
    list->Add(tree[i]);
  }

  TFile *newfile = new TFile(Form("root/ana/Sr78_ana_%s.root", type), "RECREATE");
  gROOT->cd();
  TTree *newtree = TTree::MergeTrees(list);
  newfile->cd();
  newtree->Write(0, TObject::kOverwrite);
  newfile->Close();
}

TH1 *AccTinaYY1(Int_t nBin = 90, Double_t maxTheta = 90, const char *type = "yy1", Bool_t flagDraw = false) {
  TFile *file = new TFile(Form("Analysis/Zr80_particle_%s.root", type));
  TTree *tree = static_cast<TTree *>(file->Get("PhysicsTree"));
  tree->AddFriend("SimulatedTree", Form("Simulation/Zr80_particle_%s.root", type));

  tree->Draw(Form("InitialConditions.fIC_ThetaCM>>hini(%d,0,%f)", nBin, maxTheta), 0, "GOFF");
  tree->Draw(Form("atan(sqrt(X*X+Y*Y)/Z)*TMath::RadToDeg()>>haccyy1(%d,0,%f)", nBin, maxTheta), "TrapMultiplicity > 0",
             "GOFF");

  TH1 *h1 = static_cast<TH1 *>(gDirectory->Get("haccyy1"));
  TH1 *h2 = static_cast<TH1 *>(gDirectory->Get("hini"));

  h1->Divide(h2);
  h1->SetDirectory(0);
  file->Close();
  delete file;
  if (flagDraw) h1->Draw();
  return h1;
}

TH1 *AccTinaTTT(Int_t nBin = 90, Double_t maxTheta = 90, const char *type = "yy1", Bool_t flagDraw = false) {
  TFile *file = new TFile(Form("Analysis/Zr80_particle_%s.root", type));
  TTree *tree = static_cast<TTree *>(file->Get("PhysicsTree"));
  tree->AddFriend("SimulatedTree", Form("Simulation/Zr80_particle_%s.root", type));

  tree->Draw(Form("InitialConditions.fIC_ThetaCM>>hini(%d,0,%f)", nBin, maxTheta), 0, "GOFF");
  tree->Draw(Form("atan(sqrt(X*X+Y*Y)/Z)*TMath::RadToDeg()>>haccttt(%d,0,%f)", nBin, maxTheta),
             "SquareMultiplicity > 0", "GOFF");

  TH1 *h1 = static_cast<TH1 *>(gDirectory->Get("haccttt"));
  TH1 *h2 = static_cast<TH1 *>(gDirectory->Get("hini"));

  h1->Divide(h2);
  h1->SetDirectory(0);
  file->Close();
  delete file;
  if (flagDraw) h1->Draw();
  return h1;
}

TH1 *AccTina(Int_t nBin = 90, Double_t maxTheta = 90, const char *type = "yy1") {
  TFile *file = new TFile(Form("Analysis/Zr80_particle_%s.root", type));
  TTree *tree = static_cast<TTree *>(file->Get("PhysicsTree"));
  tree->AddFriend("SimulatedTree", Form("Simulation/Zr80_particle_%s.root", type));

  tree->Draw(Form("InitialConditions.fIC_ThetaCM>>hini(%d,0,%f)", nBin, maxTheta), 0, "GOFF");
  tree->Draw(Form("atan(sqrt(X*X+Y*Y)/Z)*TMath::RadToDeg()>>hacc(%d,0,%f)", nBin, maxTheta), 0, "GOFF");

  TH1 *h1 = static_cast<TH1 *>(gDirectory->Get("hacc"));
  TH1 *h2 = static_cast<TH1 *>(gDirectory->Get("hini"));

  h1->Divide(h2);
  h1->SetDirectory(0);
  file->Close();
  delete file;

  return h1;
}

void SaveAccTina(Int_t nBin = 90, Double_t maxTheta = 90, const char *type = "yy1") {
  TH1 *hacc = AccTina(nBin, maxTheta, type);
  TH1 *haccyy1 = AccTinaYY1(nBin, maxTheta, type);
  TH1 *haccttt = AccTinaTTT(nBin, maxTheta, type);
  TFile *file = new TFile("tina_acc.root", "RECREATE");
  hacc->Write();
  haccyy1->Write();
  haccttt->Write();
  file->Close();
  delete file;
}

TH1 *LoadAccTina() {
  TFile *file = new TFile("tina_acc.root");
  TH1 *hacc = static_cast<TH1 *>(file->Get("hacc"));
  if (hacc) {
    hacc->SetDirectory(0);
    accHist = hacc;
  } else
    accHist = NULL;

  TH1 *haccyy1 = static_cast<TH1 *>(file->Get("haccyy1"));
  if (haccyy1) {
    haccyy1->SetDirectory(0);
    accYY1Hist = haccyy1;
  } else
    accYY1Hist = NULL;

  TH1 *haccttt = static_cast<TH1 *>(file->Get("haccttt"));
  if (haccttt) {
    haccttt->SetDirectory(0);
    accTTTHist = haccttt;
  } else
    accTTTHist = NULL;

  file->Close();
  delete file;
  return accHist;
}

TH1 *LoadCalcHist(const char *filename, const char *histname) {
  TFile *file = new TFile(filename);
  file->GetObject(histname, calcHist);
  if (calcHist)
    calcHist->SetDirectory(0);
  else
    calcHist = NULL;
  file->Close();
  delete file;
  return calcHist;
}

void LoadTree(const char *filename) {
  if (dataFile) delete dataFile;
  dataFile = new TFile(filename);
  dataFile->GetObject("PhysicsTree", dataTree);
}

TH1 *GetDSDW(Int_t nBin = 90, Double_t minTheta = 0, Double_t maxTheta = 90, Int_t gCut = 1, Double_t conv = 1.,
             Int_t set = 0) {
  if (!dataTree) return NULL;

  Int_t cacaoMadd;
  Double_t cacaoEdcadd[300];
  Double_t thetaL, thetaCM, eL;
  TVector3 *fragMom;
  Int_t tttMult, yy1Mult;
  fragMom = new TVector3;

  dataTree->SetBranchAddress("fragMom", &fragMom);
  dataTree->SetBranchAddress("ThetaLab", &thetaL);
  dataTree->SetBranchAddress("ThetaCM", &thetaCM);
  dataTree->SetBranchAddress("ELab", &eL);
  dataTree->SetBranchAddress("TrapMultiplicity", &yy1Mult);
  dataTree->SetBranchAddress("SquareMultiplicity", &tttMult);
  dataTree->SetBranchAddress("cacaoMadd", &cacaoMadd);
  dataTree->SetBranchAddress("cacaoEdcadd", cacaoEdcadd);

  ////////////////////////////////////////////////////////////
  // SHARAQ mom cut
  Double_t momCen = 2.03934;
  if (TString(dataFile->GetName()).Contains("45"))
    momCen = 1.86;
  else if (TString(dataFile->GetName()).Contains("40"))
    momCen = 1.74;
  else if (TString(dataFile->GetName()).Contains("55"))
    momCen = 2.11;
  else if (TString(dataFile->GetName()).Contains("60"))
    momCen = 2.22;
  ////////////////////////////////////////////////////////////

  TH1D *h1 = new TH1D(Form("hdsdw%d", set), Form("hdsdw%d", set), nBin, minTheta, maxTheta);
  TH1D *h2 = new TH1D(Form("hdsdw%d_nw", set), Form("hdsdw%d_nw", set), nBin, minTheta, maxTheta);  // no weight factor

  Int_t nEvent = 0;
  for (Int_t i = 0; i < dataTree->GetEntries(); i++) {
    dataTree->GetEntry(i);
    if (tttMult == 0 && yy1Mult == 0) continue;
    if (eL < 35) continue;                                          // ELab cut
    if (cacaoMadd < 1) continue;                                    // no gamma cut
    if (gCut == 1 && abs(cacaoEdcadd[0] - 0.471) > 0.07) continue;  // 471-keV cut
    if (gCut == 2 && abs(cacaoEdcadd[0] - 2.251) > 0.25) continue;  // 2251-keV cut
    if (gCut == 3 && abs(cacaoEdcadd[0] - 0.701) > 0.1) continue;   // 701-keV cut

    if (abs(fragMom->Mag() - momCen) > momCen * 0.03) continue;  // +-3% momentum cut

    if (set == 1 && yy1Mult > 0) continue;
    if (set == 2 && tttMult > 0) continue;

    Double_t solAngFactor = 1. / TMath::Sin(thetaL * TMath::DegToRad());
    Double_t accFactor = 1;
    //    if (accHist) {
    //      accFactor = accHist->GetBinContent(accHist->FindBin(thetaL));
    //      if (accFactor <= 0) continue;
    //      accFactor = 1./accFactor;}
    if (TString(dataFile->GetName()).Contains("tttyy1") && tttMult > 0 && thetaL < 12) continue;

    //    if (tttMult > 0) continue;
    nEvent++;
    //    cout << thetaCM << endl;
    if (tttMult > 0) {
      if (accTTTHist) {
        accFactor = accTTTHist->GetBinContent(accTTTHist->FindBin(thetaL));
        if (accFactor <= 0) continue;
        accFactor = 1. / accFactor;
      }
    } else if (yy1Mult > 0) {
      if (accYY1Hist) {
        accFactor = accYY1Hist->GetBinContent(accYY1Hist->FindBin(thetaL));
        if (accFactor <= 0) continue;
        accFactor = 1. / accFactor;
      }
    }
    h1->Fill(thetaCM, solAngFactor * accFactor);
    h2->Fill(thetaCM);
  }
  h1->Scale(1. / conv);
  std::cout << "Total event: " << nEvent << std::endl;
  return h1;
}

void DrawDSDWAndCalc(const char *filename, Int_t nBin = 90, Double_t minTheta = 0, Double_t maxTheta = 90,
                     Double_t refThetaCM = 30, Int_t gCut = 1, Double_t ymin = 1.e-4, Double_t ymax = 1,
                     Int_t set = 0) {
  LoadAccTina();

  TH1 *hc_01, *hc_21, *hc_41;
  if (TString(filename).Contains("45")) {
    hc_01 = LoadCalcHist("cs/calc_45.root", "j0_twofnr_35_43");
    hc_21 = LoadCalcHist("cs/calc_45.root", "j2_twofnr_35_43");
    hc_41 = LoadCalcHist("cs/calc_45.root", "j4_twofnr_35_43");
  }
  // hc_01 = LoadCalcHist("cs/thetacm.root","h45_0");
  //     hc_21 = LoadCalcHist("cs/thetacm.root","h45_2");
  //     hc_41 = LoadCalcHist("cs/thetacm.root","h45_4");    }
  else if (TString(filename).Contains("40")) {
    hc_01 = LoadCalcHist("cs/calc_40.root", "j0_twofnr_30_38");
    hc_21 = LoadCalcHist("cs/calc_40.root", "j2_twofnr_30_38");
    hc_41 = LoadCalcHist("cs/calc_40.root", "j4_twofnr_30_38");
  } else if (TString(filename).Contains("55")) {
    hc_01 = LoadCalcHist("cs/calc_55.root", "j0_twofnr_47_53");
    hc_21 = LoadCalcHist("cs/calc_55.root", "j2_twofnr_47_53");
    hc_41 = LoadCalcHist("cs/calc_55.root", "j4_twofnr_47_53");
  }
  //    hc_01 = LoadCalcHist("cs/thetacm.root","h55_0");
  //    hc_21 = LoadCalcHist("cs/thetacm.root","h55_2");
  //    hc_41 = LoadCalcHist("cs/thetacm.root","h55_4");    }
  else if (TString(filename).Contains("60")) {
    hc_01 = LoadCalcHist("cs/calc_60.root", "j0_twofnr_53_58");
    hc_21 = LoadCalcHist("cs/calc_60.root", "j2_twofnr_53_58");
    hc_41 = LoadCalcHist("cs/calc_60.root", "j4_twofnr_53_58");
  } else {
    hc_01 = LoadCalcHist("cs/calc.root", "j0_twofnr_42_48");
    hc_21 = LoadCalcHist("cs/calc.root", "j2_twofnr_42_48");
    hc_41 = LoadCalcHist("cs/calc.root", "j4_twofnr_42_48");
  }
  //    hc_01 = LoadCalcHist("cs/thetacm.root","h50_0");
  //    hc_21 = LoadCalcHist("cs/thetacm.root","h50_2");
  //    hc_41 = LoadCalcHist("cs/thetacm.root","h50_4");    }

  LoadTree(filename);
  TH1 *h1 = GetDSDW(nBin, minTheta, maxTheta, gCut, 1., set);
  h1->SetTitle(";#theta_{CM} [deg];# of Events (acc. & solid angle corr.)");
  h1->GetYaxis()->SetRangeUser(ymin, ymax);

  Double_t refDSDWExp = h1->GetBinContent(h1->FindBin(refThetaCM));
  Double_t refDSDWCal = hc_01->GetBinContent(hc_01->FindBin(refThetaCM));
  hc_01->Scale(refDSDWExp / refDSDWCal);
  refDSDWCal = hc_21->GetBinContent(hc_21->FindBin(refThetaCM));
  hc_21->Scale(refDSDWExp / refDSDWCal);
  refDSDWCal = hc_41->GetBinContent(hc_41->FindBin(refThetaCM));
  hc_41->Scale(refDSDWExp / refDSDWCal);

  hc_01->SetLineColor(1);
  hc_01->SetLineWidth(3);

  hc_21->SetLineColor(3);
  hc_21->SetLineWidth(3);

  hc_41->SetLineColor(4);
  hc_41->SetLineWidth(3);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetMarkerColor(2);

  h1->Draw("P, E1");
  hc_01->Draw("SAME, C, HIST");
  hc_21->Draw("SAME, C, HIST");
  hc_41->Draw("SAME, C, HIST");
  gPad->SetLogy();

  TLegend *leg = new TLegend(0.7, 0.7, 0.95, 0.95);
  leg->AddEntry(h1, "Simulation", "lp");
  leg->AddEntry(hc_01, "Calc. l = 0", "l");
  leg->AddEntry(hc_21, "Calc. l = 2", "l");
  leg->AddEntry(hc_41, "Calc. l = 4", "l");
  leg->Draw();
}

void DrawDSDWAndCalc45(const char *filename, Int_t nBin = 90, Double_t minTheta = 0, Double_t maxTheta = 90,
                       Double_t refThetaCM = 30, Int_t gCut = 1, Double_t ymin = 1.e-4, Double_t ymax = 1,
                       Int_t set = 0) {
  LoadAccTina();
  TH1 *hc_01 = LoadCalcHist("cs/calc_45.root", "j0_twofnr_35_43");
  TH1 *hc_21 = LoadCalcHist("cs/calc_45.root", "j2_twofnr_35_43");
  TH1 *hc_41 = LoadCalcHist("cs/calc_45.root", "j4_twofnr_35_43");

  LoadTree(filename);
  TH1 *h1 = GetDSDW(nBin, minTheta, maxTheta, gCut, 1., set);
  h1->SetTitle(";#theta_{CM} [deg];# of Events (acc. & solid angle corr.)");
  h1->GetYaxis()->SetRangeUser(ymin, ymax);

  Double_t refDSDWExp = h1->GetBinContent(h1->FindBin(refThetaCM));
  Double_t refDSDWCal = hc_01->GetBinContent(hc_01->FindBin(refThetaCM));
  hc_01->Scale(refDSDWExp / refDSDWCal);
  refDSDWCal = hc_21->GetBinContent(hc_21->FindBin(refThetaCM));
  hc_21->Scale(refDSDWExp / refDSDWCal);
  refDSDWCal = hc_41->GetBinContent(hc_41->FindBin(refThetaCM));
  hc_41->Scale(refDSDWExp / refDSDWCal);

  hc_01->SetLineColor(1);
  hc_01->SetLineWidth(3);

  hc_21->SetLineColor(3);
  hc_21->SetLineWidth(3);

  hc_41->SetLineColor(4);
  hc_41->SetLineWidth(3);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetMarkerColor(2);

  h1->Draw("P, E1");
  hc_01->Draw("SAME, C, HIST");
  hc_21->Draw("SAME, C, HIST");
  hc_41->Draw("SAME, C, HIST");
  gPad->SetLogy();

  TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.8);
  leg->AddEntry(h1, "Simulation", "lp");
  leg->AddEntry(hc_01, "Calc. l = 0", "l");
  leg->AddEntry(hc_21, "Calc. l = 2", "l");
  leg->AddEntry(hc_41, "Calc. l = 4", "l");
  leg->Draw();
}

void DrawGammaSpectra(const char *cut = "Max$(SquareMultiplicity > 0) == 1") {
  if (!dataTree) return;

  TH1 *h1;
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 900);
  c1->Divide(2, 2, 0.001, 0.001);

  TVirtualPad *vpad = c1->cd(1);
  dataTree->Draw("cacaoEdcadd>>hgam(250,0,1)", cut);
  h1 = static_cast<TH1 *>(gDirectory->Get("hgam"));
  h1->SetTitle(";Energy [MeV];# of Events (per 4 keV)");
  h1->GetYaxis()->SetRangeUser(0.5, 500);
  h1->SetLineColor(1);
  vpad->SetLogy();

  TBox box;
  box.SetFillStyle(3002);
  box.SetFillColor(7);
  box.DrawBox(0.239, 0.5, 0.339, 500);

  vpad = c1->cd(2);
  dataTree->Draw("cacaoEdcadd[0]:cacaoEdcadd[1]>>hgg(100,0,1,100,0,1)", cut, "COL");
  h1 = static_cast<TH1 *>(gDirectory->Get("hgg"));
  h1->SetTitle(";#gamma_{1} Energy [MeV];#gamma_{2} Energy [MeV]");

  vpad = c1->cd(3);
  dataTree->Draw("cacaoEdcadd>>hgamm(100,0,1)", cut);
  h1 = static_cast<TH1 *>(gDirectory->Get("hgamm"));
  h1->SetTitle(";Energy [MeV];# of Events (per 10 keV)");
  h1->GetYaxis()->SetRangeUser(0.5, 1500);
  h1->SetLineColor(1);
  vpad->SetLogy();

  vpad = c1->cd(4);
  TString defCut = "abs(cacaoEdcadd[1] - 0.289) < 0.05";
  if (cut) {
    defCut += " && ";
    defCut += cut;
  }
  dataTree->Draw("cacaoEdcadd[0]>>hgamgated(100,0,1)", defCut.Data());
  h1 = static_cast<TH1 *>(gDirectory->Get("hgamgated"));
  h1->SetTitle(";Energy [MeV];# of Events (per 10 keV)");
  h1->GetYaxis()->SetRangeUser(0.5, 50);
  h1->SetLineColor(1);
  vpad->SetLogy();
}

void DrawGammaSpectra2(const char *cut = "Max$(SquareMultiplicity > 0) == 1") {
  if (!dataTree) return;

  TH1 *h1;
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 900);
  c1->Divide(2, 2, 0.001, 0.001);

  TVirtualPad *vpad = c1->cd(1);
  dataTree->Draw("cacaoEdcadd>>hgam(625,0,2.5)", cut);
  h1 = static_cast<TH1 *>(gDirectory->Get("hgam"));
  h1->SetTitle(";Energy [MeV];# of Events (per 4 keV)");
  h1->SetLineColor(1);
  h1->GetYaxis()->SetRangeUser(0.5, 500);
  vpad->SetLogy();

  TBox box;
  box.SetFillStyle(3002);
  box.SetFillColor(7);
  box.DrawBox(0.239, 0.5, 0.339, 500);

  vpad = c1->cd(2);
  dataTree->Draw("cacaoEdcadd[0]:cacaoEdcadd[1]>>hgg(125,0,2.5,125,0,2.5)", cut, "COL");
  h1 = static_cast<TH1 *>(gDirectory->Get("hgg"));
  h1->SetTitle(";#gamma 1 Energy [MeV];#gamma 2 Energy [MeV]");

  vpad = c1->cd(3);
  dataTree->Draw("cacaoEdcadd>>hgamm(250,0,2.5)", cut);
  h1 = static_cast<TH1 *>(gDirectory->Get("hgamm"));
  h1->SetTitle(";Energy [MeV];# of Events (per 10 keV)");
  h1->SetLineColor(1);
  h1->GetYaxis()->SetRangeUser(0.5, 1500);
  vpad->SetLogy();

  vpad = c1->cd(4);
  TString defCut = "abs(cacaoEdcadd[1] - 0.289) < 0.05";
  if (cut) {
    defCut += " && ";
    defCut += cut;
  }
  dataTree->Draw("cacaoEdcadd[0]>>hgamgated(125,0,2.5)", defCut.Data());
  h1 = static_cast<TH1 *>(gDirectory->Get("hgamgated"));
  h1->SetTitle(";Energy [MeV];# of Events (per 20 keV)");
  h1->GetYaxis()->SetRangeUser(0.5, 50);
  h1->SetLineColor(1);
  vpad->SetLogy();
}

void DrawExThetaCM(const char *cut = NULL) {
  if (!dataTree) return;
  TText text;
  text.SetTextAlign(32);

  TH1 *h1;
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 900);
  c1->Divide(1, 2, 0.001, 0.001);

  TVirtualPad *vpad = c1->cd(1);
  dataTree->Draw("Ex>>hex(40,-1,3)", cut);
  h1 = static_cast<TH1 *>(gDirectory->Get("hex"));
  h1->SetTitle(";Ex. Spectrum [MeV]");
  h1->Fit("gaus");

  TF1 *f1 = static_cast<TF1 *>(h1->GetFunction("gaus"));
  Double_t exSigma = f1->GetParameter(2);
  Double_t exPeak = f1->GetParameter(1);
  text.DrawTextNDC(0.85, 0.8, Form("Peak = %.2f", exPeak));
  text.DrawTextNDC(0.85, 0.7, Form("FWHM = %.2f", exSigma * 2.35));

  vpad = c1->cd(2);
  dataTree->Draw("ThetaCM>>htcm(36,0,180)", cut, "E1");
  h1 = static_cast<TH1 *>(gDirectory->Get("htcm"));
  h1->SetTitle(";Theta CM [deg]");
  vpad->SetLogy();
}

void DrawExThetaCM2(const char *cut = NULL) {
  if (!dataTree) return;
  TText text;
  text.SetTextAlign(32);

  TH1 *h1;
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 900);
  c1->Divide(1, 2, 0.001, 0.001);

  TVirtualPad *vpad = c1->cd(1);
  dataTree->Draw("Ex>>hex(30,-1,5)", cut);
  h1 = static_cast<TH1 *>(gDirectory->Get("hex"));
  h1->SetTitle(";Ex. Spectrum [MeV]");
  h1->Fit("gaus");

  TF1 *f1 = h1->GetFunction("gaus");
  Double_t exSigma = f1->GetParameter(2);
  Double_t exPeak = f1->GetParameter(1);
  text.DrawTextNDC(0.85, 0.8, Form("Peak = %.2f", exPeak));
  text.DrawTextNDC(0.85, 0.7, Form("FWHM = %.2f", exSigma * 2.35));

  vpad = c1->cd(2);
  dataTree->Draw("ThetaCM>>htcm(36,0,180)", cut, "E1");
  h1 = static_cast<TH1 *>(gDirectory->Get("htcm"));
  h1->SetTitle(";Theta CM [deg]");
  vpad->SetLogy();
}

void DrawFinalBMF() {
  TCanvas *c1 = new TCanvas;
  DrawDSDWAndCalc("root/ana/Sr78_ana_cacao_ttt30cm_lh2_bmf_batch_50mev.root", 12, 6, 48, 24, 1, 0.1, 500);
  c1->Print("figs/ang_dist_bmf_cacao.pdf");
}

void DrawFinal2() {
  TCanvas *c1 = new TCanvas;
  //  DrawDSDWAndCalc("Analysis/Zr80_ana_dali_mcsm_batch.root",17,0,90,36,2,0.2,800);
  DrawDSDWAndCalc("root/ana/Sr78_ana_dali_ttt30cm_lh2_mcsm_batch_50mev.root", 12, 8, 56, 20, 2, 1, 1000);
  c1->Print("figs/ang_dist_mcsm_rev_2022.pdf");
}

void DrawGamma1() {
  TCanvas *c1 = new TCanvas;
  LoadTree("root/ana/Sr78_ana_cacao_ttt30cm_lh2_bmf_batch_50mev.root");
  DrawGammaSpectra("Max$(ELab > 35) == 1");
  gPad->GetCanvas()->Print("figs/gamma_spec_bmf_cacao.pdf");
}
void DrawGamma2() {
  TCanvas *c1 = new TCanvas;
  LoadTree("root/ana/Sr78_ana_dali_ttt30cm_lh2_mcsm_batch_50mev.root");
  DrawGammaSpectra2("Max$(ELab > 35) == 1");
  gPad->GetCanvas()->Print("figs/gamma_spec_mcsm_cacao.pdf");
}

void DrawMD(const char *filename, Int_t nBin = 90, Double_t minTheta = 0, Double_t maxTheta = 90,
            Double_t refThetaCM = 30, Int_t gCut = 1, Double_t ymin = 1.e-4, Double_t ymax = 1) {
  LoadAccTina();

  TH1 *hc_01, *hc_21, *hc_41;
  if (TString(filename).Contains("45")) {
    hc_01 = LoadCalcHist("cs/calc_45.root", "j0_twofnr_35_43");
    hc_21 = LoadCalcHist("cs/calc_45.root", "j2_twofnr_35_43");
    hc_41 = LoadCalcHist("cs/calc_45.root", "j4_twofnr_35_43");
  } else if (TString(filename).Contains("40")) {
    hc_01 = LoadCalcHist("cs/calc_40.root", "j0_twofnr_30_38");
    hc_21 = LoadCalcHist("cs/calc_40.root", "j2_twofnr_30_38");
    hc_41 = LoadCalcHist("cs/calc_40.root", "j4_twofnr_30_38");
  } else if (TString(filename).Contains("55")) {
    hc_01 = LoadCalcHist("cs/calc_55.root", "j0_twofnr_47_53");
    hc_21 = LoadCalcHist("cs/calc_55.root", "j2_twofnr_47_53");
    hc_41 = LoadCalcHist("cs/calc_55.root", "j4_twofnr_47_53");
  } else if (TString(filename).Contains("60")) {
    hc_01 = LoadCalcHist("cs/calc_60.root", "j0_twofnr_53_58");
    hc_21 = LoadCalcHist("cs/calc_60.root", "j2_twofnr_53_58");
    hc_41 = LoadCalcHist("cs/calc_60.root", "j4_twofnr_53_58");
  } else {
    hc_01 = LoadCalcHist("cs/calc.root", "j0_twofnr_42_48");
    hc_21 = LoadCalcHist("cs/calc.root", "j2_twofnr_42_48");
    hc_41 = LoadCalcHist("cs/calc.root", "j4_twofnr_42_48");
  }

  LoadTree(filename);
  TH1 *h1 = GetDSDW(nBin, minTheta, maxTheta, gCut);
  h1->SetTitle(";#theta_{CM} [deg];# of Events (acc. & solid angle corr.)");
  h1->GetYaxis()->SetRangeUser(ymin, ymax);

  Double_t refScale = 100;
  Double_t portion_bmf[3] = {0.89, 0.1, 0.01};
  Double_t portion_mcsm[3] = {0.62, 0.32, 0.06};
  Double_t *portion;
  if (TString(filename).Contains("bmf"))
    portion = portion_bmf;
  else
    portion = portion_mcsm;

  hc_01->Scale(refScale * portion[0] / hc_01->Integral());
  hc_21->Scale(refScale * portion[1] / hc_21->Integral());
  hc_41->Scale(refScale * portion[2] / hc_41->Integral());

  TH1 *hc = static_cast<TH1 *>(hc_01->Clone("h_md"));
  hc->Add(hc_21);
  hc->Add(hc_41);

  Double_t refDSDWExp = h1->GetBinContent(h1->FindBin(refThetaCM));
  Double_t refDSDWCal = hc->GetBinContent(hc->FindBin(refThetaCM));
  hc->Scale(refDSDWExp / refDSDWCal);
  hc->SetLineColor(2);
  hc->SetLineWidth(3);

  hc_01->Scale(refDSDWExp / refDSDWCal);
  hc_21->Scale(refDSDWExp / refDSDWCal);
  hc_41->Scale(refDSDWExp / refDSDWCal);

  hc_01->SetLineColor(1);
  hc_01->SetLineStyle(2);
  hc_01->SetLineWidth(3);

  hc_21->SetLineColor(3);
  hc_21->SetLineStyle(3);
  hc_21->SetLineWidth(3);

  hc_41->SetLineColor(4);
  hc_41->SetLineStyle(4);
  hc_41->SetLineWidth(3);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetMarkerColor(2);

  h1->Draw("P, E1");
  hc->Draw("SAME, C, HIST");
  hc_01->Draw("SAME, C, HIST");
  hc_21->Draw("SAME, C, HIST");
  hc_41->Draw("SAME, C, HIST");
  gPad->SetLogy();

  TLegend *leg = new TLegend(0.7, 0.7, 0.95, 0.95);
  leg->AddEntry(h1, "Simulation", "lp");
  leg->AddEntry(hc, "Mult.Decomp.", "l");
  leg->AddEntry(hc_01, "l = 0", "l");
  leg->AddEntry(hc_21, "l = 2", "l");
  leg->AddEntry(hc_41, "l = 4", "l");
  leg->Draw();
}

void CheckSHARAQAcc(const char *filename) {
  LoadTree(filename);
  dataTree->Draw("fragMom.Mag()>>hf(1000,1.3,3.3)");
  TH1 *h1 = static_cast<TH1 *>(gDirectory->Get("hf"));
  h1->Fit("gaus");
  Double_t cen = h1->GetFunction("gaus")->GetParameter(1);
  Int_t cenBin = h1->FindBin(cen);
  Int_t uppBin = h1->FindBin(cen * 1.03);
  Int_t lowBin = h1->FindBin(cen * 0.97);
  Int_t accEntries = h1->Integral(lowBin, uppBin);
  Int_t totEntries = h1->Integral();
  vector<double> acc;
  for (Int_t i = -100; i <= 100; i++) {
    acc.push_back((double)h1->Integral(lowBin + i, uppBin + i) / h1->Integral());
  }
  Int_t maxInd = max_element(acc.begin(), acc.end()) - acc.begin();
  Double_t center = h1->GetBinCenter(cenBin + maxInd - 100);
  Double_t lowLim = h1->GetBinCenter(lowBin + maxInd - 100);
  Double_t uppLim = h1->GetBinCenter(uppBin + maxInd - 100);

  gPad->GetCanvas()->Update();
  TLine line;
  line.DrawLine(lowLim, 0, lowLim, gPad->GetUymax());
  line.DrawLine(uppLim, 0, uppLim, gPad->GetUymax());
  std::cout << center << ", " << acc[maxInd] << std::endl;
}

void DrawFinalMD1() {
  TCanvas *c1 = new TCanvas;
  //  DrawDSDWAndCalc("Analysis/Zr80_ana_dali_bmf_batch.root",18,0,90,37,1,0.2,800);
  // DrawDSDWAndCalc("Analysis/Zr80_ana_dali_ttt30cm_lh2_bmf_batch_50mev.root",12,8,56,23,1,1,1000);
  // DrawMD("root/ana/Sr78_ana_cacao_ttt30cm_lh2_bmf_batch_50mev.root", 12, 6, 48, 24, 1, 0.1, 500);
  DrawMD("root/ana/Sr78_ana_cacao_ttt30cm_lh2_bmf_50mev.root", 12, 6, 48, 24, 1, 0.1, 500);
  c1->Print("figs/ang_dist_bmf_cacao_md.pdf");
}

void DrawFinalMD2() {
  TCanvas *c1 = new TCanvas;
  //  DrawDSDWAndCalc("Analysis/Zr80_ana_dali_mcsm_batch.root",17,0,90,36,2,0.2,800);
  //  DrawDSDWAndCalc("Analysis/Zr80_ana_dali_ttt30cm_lh2_mcsm_batch_50mev.root",12,8,56,20,2,1,1000);
  DrawMD("root/ana/Sr78_ana_cacao_ttt30cm_lh2_mcsm_50mev.root", 12, 6, 48, 30, 2, 0.1, 500);
  c1->Print("figs/ang_dist_mcsm_cacao_md.pdf");
}
