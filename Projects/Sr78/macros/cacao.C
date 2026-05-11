#ifndef __CLING__
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TText.h"
#include "TTree.h"
#endif

using namespace std;

void cacao() {
    bool flagTest = false;
    int nBins = 750;

    vector<int> ens = {471, 2251};
    vector<int> dists = {105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205};
    vector<vector<Double_t>> effs;
    for (const auto &dist : dists) {
        vector<Double_t> ns;
        for (const auto &en : ens) {
            TFile *file = new TFile(Form("root/ana/cacao_ana_%d_%d.root", en, dist));
            TTree *tree = static_cast<TTree *>(file->Get("PhysicsTree"));
            Int_t nEntries = tree->GetEntries();
            TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
            TH1F *h1 = new TH1F("h1", "h1", nBins, 0, 3);
            tree->Draw("cacaoEadd>>h1");
            Double_t low_lim = en < 1000 ? 0.4 : 2.1;
            Double_t upp_lim = en < 1000 ? 0.6 : 2.4;
            h1->Fit("gausn", 0, 0, low_lim, upp_lim);
            TF1 *f1 = h1->GetFunction("gausn");

            Double_t pars[6];
            for (int i = 0; i < 3; i++) pars[i] = f1->GetParameter(i);
            for (int i = 3; i < 6; i++) pars[i] = 0;

            TF1 *f2 = new TF1("f2", "gausn(0)+pol2(3)", 0, 3);
            f2->SetParameters(pars);
            h1->Fit(f2, 0, 0, low_lim, upp_lim);
            Double_t n = f2->GetParameter(0) / h1->GetBinWidth(1);
            Double_t eff = n / nEntries;
            ns.push_back(eff);
            c1->SaveAs(Form("figs/cacao/cacao_and_%d_%d.pdf", en, dist));
            delete c1;
            file->Close();
            delete file;
            if (flagTest) break;
        }
        effs.push_back(ns);
        if (flagTest) break;
    }
    for (const auto &ns : effs) {
        for (const auto &n : ns) {
            cout << n << " ";
        }
        cout << endl;
    }
    // effs 벡터를 이용하여 distance에 따른 efficiency 변화를 TGraph로 그립니다.
    for (size_t i = 0; i < ens.size(); ++i) {
        std::vector<Double_t> y;
        for (size_t j = 0; j < dists.size(); ++j) {
            y.push_back(effs[j][i]);
        }
        TGraph *gr = new TGraph(dists.size());
        for (size_t j = 0; j < dists.size(); ++j) {
            gr->SetPoint(j, dists[j], y[j]);
        }
        gr->SetTitle(Form("Efficiency vs Distance (E = %d keV);Distance (mm);Efficiency", ens[i]));
        gr->SetMarkerStyle(20 + i);
        gr->SetMarkerColor(2 + i);
        gr->SetLineColor(2 + i);
        TCanvas *c2 = new TCanvas(Form("c2_%d", ens[i]), Form("Efficiency_%d", ens[i]), 800, 600);
        gr->Draw("APL");
        // 105 mm에 해당하는 포인트 옆에 "Default"라는 텍스트를 추가
        for (size_t j = 0; j < dists.size(); ++j) {
            if (dists[j] == 105) {
                double x = dists[j] + 2;
                double yval = y[j];
                TLatex *latex = new TLatex(x, yval, "Default");
                latex->SetTextSize(0.03);
                latex->SetTextAlign(12);  // 왼쪽 정렬
                latex->Draw();
            }
        }
        c2->SaveAs(Form("figs/cacao/eff_vs_dist_%d.pdf", ens[i]));
        delete c2;
        delete gr;
    }
}