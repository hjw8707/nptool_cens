#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TText.h"
#include "TTree.h"
#include "TVirtualPad.h"

using namespace std;

void DrawEachAll(const char* type, Double_t range = 0.06);
void DrawEachDG(const char* type, const char* nn, Double_t fitPos, Double_t range = 0.06);
void DrawEachDGMCSM(const char* type, const char* nn, Double_t fitPos, Double_t range = 0.06);
void DrawEach(const char* type, const char* file, Double_t fitPos, Double_t range = 0.06);
void DrawEachMCSM(const char* type, const char* file, Double_t fitPos, Double_t range = 0.06);
void LoadAll();
void DrawCoin();
void DrawEachAddback(const char* type, const char* file);
void DrawEachAddbackAll(const char* type);

const char* cut = "ELab > 35";
TFile* file = NULL;
TTree* tree = NULL;

TDirectory* tdir = NULL;

void LoadAll() {
  tdir = new TDirectory("tdir", "tdir");

  const char* files[] = {"1st0", "2nd0", "1st2", "2nd2"};
  for (int i = 0; i < 4; i++) {
    TH1* h1;

    if (file) {
      file->Close();
      delete file;
    }
    file = new TFile(Form("root/ana/Sr78_ana_%s.root", files[i]));
    tree = static_cast<TTree*>(file->Get("PhysicsTree"));

    tree->Draw("Ex>>hex(120,-2,4)", cut);
    h1 = static_cast<TH1*>(gDirectory->Get("hex"));
    h1->SetTitle(";Ex [MeV];Event #");
    h1 = static_cast<TH1*>(h1->Clone(Form("hex_%s", files[i])));
    h1->SetDirectory(tdir);

    tree->Draw("ThetaCM>>htcm(18,0,90)", cut, "E1");
    h1 = static_cast<TH1*>(gDirectory->Get("htcm"));
    h1->SetTitle(";Theta CM [deg];Event #");
    h1 = static_cast<TH1*>(h1->Clone(Form("htcm_%s", files[i])));
    h1->SetDirectory(tdir);

    tree->Draw("Edc_add>>hgam(100,0,1)", cut);
    h1 = static_cast<TH1*>(gDirectory->Get("hgam"));
    h1->SetTitle(";D.C. Gamma [MeV];Event #");
    h1 = static_cast<TH1*>(h1->Clone(Form("hgam_%s", files[i])));
    h1->SetDirectory(tdir);
  }

  tdir->cd();
}

void DrawEachAll(const char* type, Double_t range) {
  const char* files[] = {"1st0", "1st2", "2nd0", "2nd2"};
  Double_t fPos[4] = {0., 0.289, 0.471, 0.811};
  for (int i = 0; i < 4; i++) {
    DrawEach(type, files[i], fPos[i], range);
    gPad->GetCanvas()->Update();
    getchar();
  }
}

void DrawEachAllMCSM(const char* type, Double_t range) {
  const char* files[] = {"1st0", "1st2", "2nd0", "2nd2"};
  Double_t fPos[4] = {0., 0.289, 2.251, 0.701};
  for (int i = 0; i < 4; i++) {
    DrawEachMCSM(type, files[i], fPos[i], range);
    gPad->GetCanvas()->Update();
    getchar();
  }
}

void DrawEachAddbackAll(const char* type) {
  const char* files[] = {"1st0", "1st2", "2nd0", "2nd2"};
  for (int i = 0; i < 4; i++) {
    DrawEachAddback(type, files[i]);
    gPad->GetCanvas()->Update();
    getchar();
  }
}

void DrawEachAllDG(const char* type, Double_t range) {
  const char* files[] = {"1st0", "1st2", "2nd0", "2nd2"};
  Double_t fPos[4] = {0., 0.289, 0.471, 0.811};
  for (int i = 0; i < 4; i++) {
    DrawEachDG(type, files[i], fPos[i], range);
    gPad->GetCanvas()->Update();
    getchar();
  }
}

void DrawEachAllDGMCSM(const char* type, Double_t range) {
  const char* files[] = {"1st0", "1st2", "2nd0", "2nd2"};
  Double_t fPos[4] = {0., 0.289, 2.251, 0.701};
  for (int i = 0; i < 4; i++) {
    DrawEachDGMCSM(type, files[i], fPos[i], range);
    gPad->GetCanvas()->Update();
    getchar();
  }
}

void DrawEach(const char* type, const char* nn, Double_t fitPos, Double_t range) {
  TText text;
  text.SetTextAlign(32);
  if (file) {
    file->Close();
    delete file;
  }

  file = new TFile(Form("root/ana/Sr78_ana_%s_%s.root", nn, type));
  tree = static_cast<TTree*>(file->Get("PhysicsTree"));

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
  c1->Divide(2, 2);

  TH1* h1;
  TF1* f1;

  c1->cd(1);
  tree->Draw("Ex>>hex(120,-2,4)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hex"));
  h1->SetTitle(";Ex [MeV];Event #");
  h1->Fit("gaus");
  f1 = h1->GetFunction("gaus");
  Double_t exSigma = f1->GetParameter(2);
  Double_t nEventEx = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventEx));
  text.DrawTextNDC(0.85, 0.7, Form("FWHM = %.2f", exSigma * 2.35));

  TVirtualPad* vpad = c1->cd(2);
  tree->Draw("ThetaCM>>htcm(18,0,90)", cut, "E1");
  h1 = static_cast<TH1*>(gDirectory->Get("htcm"));
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetTitle(";Theta CM [deg];Event #");
  Double_t nEventTh = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventTh));
  vpad->SetLogy();

  Double_t fitRange = range;
  Double_t fitCor = 1.015;
  TString gTypeC(type);
  TString gType;
  if (gTypeC.BeginsWith("dali"))
    gType = "dali";
  else if (gTypeC.BeginsWith("cacao"))
    gType = "cacao";
  else if (gTypeC.BeginsWith("csi"))
    gType = "csi";
  else
    gType = "grape";
  c1->cd(3);
  tree->Draw(Form("%sEdcadd>>hgam(1000,0,1)", gType.Data()));
  h1 = static_cast<TH1*>(gDirectory->Get("hgam"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  h1->Fit("gausn", 0, 0, fitPos * (fitCor - fitRange), fitPos * (fitCor + fitRange));
  if (h1->GetFunction("gausn")) {
    f1 = h1->GetFunction("gausn");
    Double_t nPeak = f1->GetParameter(0) / h1->GetBinWidth(1);
    text.DrawTextNDC(0.85, 0.8, Form("Peak Ev. = %.0f", nPeak));
  }

  c1->cd(4);
  tree->Draw(Form("%sEdcadd>>hgamcut(1000,0,1)", gType.Data()), cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  h1->Fit("gausn", 0, 0, fitPos * (fitCor - fitRange), fitPos * (fitCor + fitRange));
  if (h1->GetFunction("gausn")) {
    f1 = h1->GetFunction("gausn");
    Double_t nPeak = f1->GetParameter(0) / h1->GetBinWidth(1);
    text.DrawTextNDC(0.85, 0.8, Form("Peak Ev. = %.0f", nPeak));
  }

  c1->Update();
  c1->Print(Form("figs/sum_%s_%s.pdf", nn, type));
}

void DrawEachMCSM(const char* type, const char* nn, Double_t fitPos, Double_t range) {
  TText text;
  text.SetTextAlign(32);
  if (file) {
    file->Close();
    delete file;
  }

  file = new TFile(Form("root/ana/Sr78_ana_%s_%s.root", nn, type));
  tree = static_cast<TTree*>(file->Get("PhysicsTree"));

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
  c1->Divide(2, 2);

  TH1* h1;
  TF1* f1;

  c1->cd(1);
  tree->Draw("Ex>>hex(160,-2,6)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hex"));
  h1->SetTitle(";Ex [MeV];Event #");
  h1->Fit("gaus");
  f1 = h1->GetFunction("gaus");
  Double_t exSigma = f1->GetParameter(2);
  Double_t nEventEx = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventEx));
  text.DrawTextNDC(0.85, 0.7, Form("FWHM = %.2f", exSigma * 2.35));

  TVirtualPad* vpad = c1->cd(2);
  tree->Draw("ThetaCM>>htcm(18,0,90)", cut, "E1");
  h1 = static_cast<TH1*>(gDirectory->Get("htcm"));
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetTitle(";Theta CM [deg];Event #");
  Double_t nEventTh = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventTh));
  vpad->SetLogy();

  Double_t fitRange = range;
  Double_t fitCor = 1.015;
  TString gTypeC(type);
  TString gType;
  if (gTypeC.BeginsWith("dali"))
    gType = "dali";
  else if (gTypeC.BeginsWith("cacao"))
    gType = "cacao";
  else
    gType = "grape";
  c1->cd(3);
  tree->Draw(Form("%sEdcadd>>hgam(2500,0,2.5)", gType.Data()));
  h1 = static_cast<TH1*>(gDirectory->Get("hgam"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  h1->Fit("gausn", 0, 0, fitPos * (fitCor - fitRange), fitPos * (fitCor + fitRange));
  if (h1->GetFunction("gausn")) {
    f1 = h1->GetFunction("gausn");
    Double_t nPeak = f1->GetParameter(0) / h1->GetBinWidth(1);
    text.DrawTextNDC(0.85, 0.8, Form("Peak Ev. = %.0f", nPeak));
  }

  c1->cd(4);
  tree->Draw(Form("%sEdcadd>>hgamcut(2500,0,2.5)", gType.Data()), cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  h1->Fit("gausn", 0, 0, fitPos * (fitCor - fitRange), fitPos * (fitCor + fitRange));
  if (h1->GetFunction("gausn")) {
    f1 = h1->GetFunction("gausn");
    Double_t nPeak = f1->GetParameter(0) / h1->GetBinWidth(1);
    text.DrawTextNDC(0.85, 0.8, Form("Peak Ev. = %.0f", nPeak));
  }

  c1->Update();
  c1->Print(Form("figs/sum_%s_%s.pdf", nn, type));
}

void DrawEachDG(const char* type, const char* nn, Double_t fitPos, Double_t range) {
  TText text;
  text.SetTextAlign(32);
  if (file) {
    file->Close();
    delete file;
  }

  file = new TFile(Form("root/ana/Sr78_ana_%s_%s.root", nn, type));
  tree = static_cast<TTree*>(file->Get("PhysicsTree"));

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
  c1->Divide(2, 2);

  TH1* h1;
  TF1* f1;

  c1->cd(1);
  tree->Draw("Ex>>hex(120,-2,4)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hex"));
  h1->SetTitle(";Ex [MeV];Event #");
  h1->Fit("gaus");
  f1 = h1->GetFunction("gaus");
  Double_t exSigma = f1->GetParameter(2);
  Double_t nEventEx = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventEx));
  text.DrawTextNDC(0.85, 0.7, Form("FWHM = %.2f", exSigma * 2.35));

  TVirtualPad* vpad = c1->cd(2);
  tree->Draw("ThetaCM>>htcm(18,0,90)", cut, "E1");
  h1 = static_cast<TH1*>(gDirectory->Get("htcm"));
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetTitle(";Theta CM [deg];Event #");
  Double_t nEventTh = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventTh));
  vpad->SetLogy();

  Double_t fitRange = range;
  Double_t fitCor = 1.015;
  c1->cd(3);
  tree->Draw("grapeEdcadd>>hgamg(1200,0,1.2)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamg"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  h1->Fit("gausn", 0, 0, fitPos * (fitCor - fitRange), fitPos * (fitCor + fitRange));
  if (h1->GetFunction("gausn")) {
    f1 = h1->GetFunction("gausn");
    Double_t nPeak = f1->GetParameter(0) / h1->GetBinWidth(1);
    text.DrawTextNDC(0.85, 0.8, Form("Peak Ev. = %.0f", nPeak));
  }

  c1->cd(4);
  tree->Draw("daliEdcadd>>hgamd(1200,0,1.2)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamd"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  h1->Fit("gausn", 0, 0, fitPos * (fitCor - fitRange), fitPos * (fitCor + fitRange));
  if (h1->GetFunction("gausn")) {
    f1 = h1->GetFunction("gausn");
    Double_t nPeak = f1->GetParameter(0) / h1->GetBinWidth(1);
    text.DrawTextNDC(0.85, 0.8, Form("Peak Ev. = %.0f", nPeak));
  }

  c1->Update();
  c1->Print(Form("figs/sum_%s_%s.pdf", nn, type));
}

void DrawEachDGMCSM(const char* type, const char* nn, Double_t fitPos, Double_t range) {
  TText text;
  text.SetTextAlign(32);
  if (file) {
    file->Close();
    delete file;
  }

  file = new TFile(Form("root/ana/Sr78_ana_%s_%s.root", nn, type));
  tree = static_cast<TTree*>(file->Get("PhysicsTree"));

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
  c1->Divide(2, 2);

  TH1* h1;
  TF1* f1;

  c1->cd(1);
  tree->Draw("Ex>>hex(120,-2,4)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hex"));
  h1->SetTitle(";Ex [MeV];Event #");
  h1->Fit("gaus");
  f1 = h1->GetFunction("gaus");
  Double_t exSigma = f1->GetParameter(2);
  Double_t nEventEx = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventEx));
  text.DrawTextNDC(0.85, 0.7, Form("FWHM = %.2f", exSigma * 2.35));

  TVirtualPad* vpad = c1->cd(2);
  tree->Draw("ThetaCM>>htcm(18,0,90)", cut, "E1");
  h1 = static_cast<TH1*>(gDirectory->Get("htcm"));
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetTitle(";Theta CM [deg];Event #");
  Double_t nEventTh = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventTh));
  vpad->SetLogy();

  Double_t fitRange = range;
  Double_t fitCor = 1.015;
  c1->cd(3);
  tree->Draw("grapeEdcadd>>hgamg(2500,0,2.5)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamg"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  h1->Fit("gausn", 0, 0, fitPos * (fitCor - fitRange), fitPos * (fitCor + fitRange));
  if (h1->GetFunction("gausn")) {
    f1 = h1->GetFunction("gausn");
    Double_t nPeak = f1->GetParameter(0) / h1->GetBinWidth(1);
    text.DrawTextNDC(0.85, 0.8, Form("Peak Ev. = %.0f", nPeak));
  }

  c1->cd(4);
  tree->Draw("daliEdcadd>>hgamd(2500,0,2.5)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamd"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  h1->Fit("gausn", 0, 0, fitPos * (fitCor - fitRange), fitPos * (fitCor + fitRange));
  if (h1->GetFunction("gausn")) {
    f1 = h1->GetFunction("gausn");
    Double_t nPeak = f1->GetParameter(0) / h1->GetBinWidth(1);
    text.DrawTextNDC(0.85, 0.8, Form("Peak Ev. = %.0f", nPeak));
  }

  c1->Update();
  c1->Print(Form("figs/sum_%s_%s.pdf", nn, type));
}

void DrawCoin() {
  TText text;
  text.SetTextAlign(32);
  if (file) {
    file->Close();
    delete file;
  }

  file = new TFile("root/ana/Sr78_ana.root");
  tree = static_cast<TTree*>(file->Get("PhysicsTree"));

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
  c1->Divide(2, 2);

  TH1* h1;

  c1->cd(1);
  tree->Draw("Ex>>hex(120,-2,4)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hex"));
  h1->SetTitle(";Ex [MeV];Event #");
  //  h1->Fit("gaus");
  Double_t nEventEx = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventEx));

  tree->Draw("Ex>>hex1(120,-2,4)", cut, "same", 100000, 0);
  h1 = static_cast<TH1*>(gDirectory->Get("hex1"));
  h1->SetLineColor(2);
  h1->SetLineStyle(2);

  tree->Draw("Ex>>hex2(120,-2,4)", cut, "same", 100000, 100000);
  h1 = static_cast<TH1*>(gDirectory->Get("hex2"));
  h1->SetLineColor(3);
  h1->SetLineStyle(2);

  tree->Draw("Ex>>hex3(120,-2,4)", cut, "same", 100000, 200000);
  h1 = static_cast<TH1*>(gDirectory->Get("hex3"));
  h1->SetLineColor(4);
  h1->SetLineStyle(2);

  tree->Draw("Ex>>hex4(120,-2,4)", cut, "same", 100000, 300000);
  h1 = static_cast<TH1*>(gDirectory->Get("hex4"));
  h1->SetLineColor(6);
  h1->SetLineStyle(2);

  c1->cd(2);
  tree->Draw("ThetaCM>>htcm(18,0,90)", cut, "E1");
  h1 = static_cast<TH1*>(gDirectory->Get("htcm"));
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetTitle(";Theta CM [deg];Event #");
  Double_t nEventTh = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventTh));

  tree->Draw("ThetaCM>>htcm1(18,0,90)", cut, "same", 100000, 0);
  h1 = static_cast<TH1*>(gDirectory->Get("htcm1"));
  h1->SetLineColor(2);
  h1->SetLineStyle(2);

  tree->Draw("ThetaCM>>htcm2(18,0,90)", cut, "same", 100000, 100000);
  h1 = static_cast<TH1*>(gDirectory->Get("htcm2"));
  h1->SetLineColor(3);
  h1->SetLineStyle(2);

  tree->Draw("ThetaCM>>htcm3(18,0,90)", cut, "same", 100000, 200000);
  h1 = static_cast<TH1*>(gDirectory->Get("htcm3"));
  h1->SetLineColor(4);
  h1->SetLineStyle(2);

  tree->Draw("ThetaCM>>htcm4(18,0,90)", cut, "same", 100000, 300000);
  h1 = static_cast<TH1*>(gDirectory->Get("htcm4"));
  h1->SetLineColor(6);
  h1->SetLineStyle(2);

  c1->cd(3);
  tree->Draw("Edc_add>>hgam(100,0,1)");
  h1 = static_cast<TH1*>(gDirectory->Get("hgam"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");

  tree->Draw("Edc_add>>hgam1(100,0,1)", 0, "same", 100000, 0);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam1"));
  h1->SetLineColor(2);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgam2(100,0,1)", 0, "same", 100000, 100000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam2"));
  h1->SetLineColor(3);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgam3(100,0,1)", 0, "same", 100000, 200000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam3"));
  h1->SetLineColor(4);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgam4(100,0,1)", 0, "same", 100000, 300000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam4"));
  h1->SetLineColor(6);
  h1->SetLineStyle(2);

  c1->cd(4);
  tree->Draw("Edc_add>>hgamcut(100,0,1)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");

  tree->Draw("Edc_add>>hgamcut1(100,0,1)", cut, "same", 100000, 0);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut1"));
  h1->SetLineColor(2);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgamcut2(100,0,1)", cut, "same", 100000, 100000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut2"));
  h1->SetLineColor(3);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgamcut3(100,0,1)", cut, "same", 100000, 200000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut3"));
  h1->SetLineColor(4);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgamcut4(100,0,1)", cut, "same", 100000, 300000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut4"));
  h1->SetLineColor(6);
  h1->SetLineStyle(2);

  c1->Update();
  c1->Print("figs/sum_coin.pdf");
}

void DrawCoinBatch() {
  TText text;
  text.SetTextAlign(32);
  if (file) {
    file->Close();
    delete file;
  }

  file = new TFile("root/ana/Sr78_ana_batch.root");
  tree = static_cast<TTree*>(file->Get("PhysicsTree"));

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
  c1->Divide(2, 2);

  TH1* h1;

  c1->cd(1);
  tree->Draw("Ex>>hex(120,-2,4)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hex"));
  h1->SetTitle(";Ex [MeV];Event #");
  //  h1->Fit("gaus");
  Double_t nEventEx = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventEx));

  tree->Draw("Ex>>hex1(120,-2,4)", cut, "same", 440000, 0);
  h1 = static_cast<TH1*>(gDirectory->Get("hex1"));
  h1->SetLineColor(2);
  h1->SetLineStyle(2);

  tree->Draw("Ex>>hex2(120,-2,4)", cut, "same", 87000, 440000);
  h1 = static_cast<TH1*>(gDirectory->Get("hex2"));
  h1->SetLineColor(3);
  h1->SetLineStyle(2);

  tree->Draw("Ex>>hex3(120,-2,4)", cut, "same", 12000, 527000);
  h1 = static_cast<TH1*>(gDirectory->Get("hex3"));
  h1->SetLineColor(4);
  h1->SetLineStyle(2);

  tree->Draw("Ex>>hex4(120,-2,4)", cut, "same", 7000, 539000);
  h1 = static_cast<TH1*>(gDirectory->Get("hex4"));
  h1->SetLineColor(6);
  h1->SetLineStyle(2);

  c1->cd(2);
  tree->Draw("ThetaCM>>htcm(18,0,90)", cut, "E1");
  h1 = static_cast<TH1*>(gDirectory->Get("htcm"));
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.0);
  h1->SetTitle(";Theta CM [deg];Event #");
  Double_t nEventTh = h1->Integral();
  text.DrawTextNDC(0.85, 0.8, Form("Event = %.0f", nEventTh));

  tree->Draw("ThetaCM>>htcm1(18,0,90)", cut, "same", 440000, 0);
  h1 = static_cast<TH1*>(gDirectory->Get("htcm1"));
  h1->SetLineColor(2);
  h1->SetLineStyle(2);

  tree->Draw("ThetaCM>>htcm2(18,0,90)", cut, "same", 87000, 440000);
  h1 = static_cast<TH1*>(gDirectory->Get("htcm2"));
  h1->SetLineColor(3);
  h1->SetLineStyle(2);

  tree->Draw("ThetaCM>>htcm3(18,0,90)", cut, "same", 12000, 527000);
  h1 = static_cast<TH1*>(gDirectory->Get("htcm3"));
  h1->SetLineColor(4);
  h1->SetLineStyle(2);

  tree->Draw("ThetaCM>>htcm4(18,0,90)", cut, "same", 7000, 539000);
  h1 = static_cast<TH1*>(gDirectory->Get("htcm4"));
  h1->SetLineColor(6);
  h1->SetLineStyle(2);

  c1->cd(3);
  tree->Draw("Edc_add>>hgam(100,0,1)");
  h1 = static_cast<TH1*>(gDirectory->Get("hgam"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");

  tree->Draw("Edc_add>>hgam1(100,0,1)", 0, "same", 440000, 0);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam1"));
  h1->SetLineColor(2);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgam2(100,0,1)", 0, "same", 87000, 440000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam2"));
  h1->SetLineColor(3);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgam3(100,0,1)", 0, "same", 12000, 527000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam3"));
  h1->SetLineColor(4);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgam4(100,0,1)", 0, "same", 7000, 539000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam4"));
  h1->SetLineColor(6);
  h1->SetLineStyle(2);

  c1->cd(4);
  tree->Draw("Edc_add>>hgamcut(100,0,1)", cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");

  tree->Draw("Edc_add>>hgamcut1(100,0,1)", cut, "same", 440000, 0);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut1"));
  h1->SetLineColor(2);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgamcut2(100,0,1)", cut, "same", 87000, 440000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut2"));
  h1->SetLineColor(3);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgamcut3(100,0,1)", cut, "same", 12000, 527000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut3"));
  h1->SetLineColor(4);
  h1->SetLineStyle(2);

  tree->Draw("Edc_add>>hgamcut4(100,0,1)", cut, "same", 7000, 539000);
  h1 = static_cast<TH1*>(gDirectory->Get("hgamcut4"));
  h1->SetLineColor(6);
  h1->SetLineStyle(2);

  c1->Update();
  c1->Print("figs/sum_coin_batch.pdf");
}

void DrawEachAddback(const char* type, const char* nn) {
  TText text;
  text.SetTextAlign(32);
  if (file) {
    file->Close();
    delete file;
  }

  file = new TFile(Form("root/ana/Sr78_ana_%s_%s.root", nn, type));
  tree = static_cast<TTree*>(file->Get("PhysicsTree"));

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
  c1->Divide(2, 2);

  TH1* h1;
  TF1* f1;

  c1->cd(1);
  tree->Draw(Form("%sEdcadd[0]>>hgam0(1200,0,1.2)", type), cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam0"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  text.DrawTextNDC(0.85, 0.8, "Highest E");

  c1->cd(2);
  tree->Draw(Form("%sEdcadd[1]>>hgam1(1200,0,1.2)", type), cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam1"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  text.DrawTextNDC(0.85, 0.8, "2nd Highest E");

  c1->cd(3);
  tree->Draw(Form("%sEdcadd[2]>>hgam2(1200,0,1.2)", type), cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam2"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  text.DrawTextNDC(0.85, 0.8, "2nd Highest E");

  c1->cd(4);
  tree->Draw(Form("%sEdcadd[3]>>hgam3(1200,0,1.2)", type), cut);
  h1 = static_cast<TH1*>(gDirectory->Get("hgam3"));
  h1->SetTitle(";D.C. Gamma [MeV];Event #");
  text.DrawTextNDC(0.85, 0.8, "2nd Highest E");

  c1->Update();
  c1->Print(Form("figs/sumab_%s_%s.pdf", nn, type));
}
