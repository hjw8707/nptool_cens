void plot_hpge_scan_overlay(const char* inputDir = "root/ana",
                            const char* prefix = "background_plus_co60_x",
                            const char* outputBase = "root/ana/background_plus_co60_overlay") {
  gROOT->SetBatch(kTRUE);

  TSystemDirectory dir("hpge_ana", inputDir);
  TList* files = dir.GetListOfFiles();
  if (!files) {
    Error("plot_hpge_scan_overlay", "Cannot list %s", inputDir);
    return;
  }

  std::vector<TString> names;
  TIter next(files);
  while (TObject* obj = next()) {
    TString name = obj->GetName();
    if (name.BeginsWith(prefix) && name.EndsWith(".root"))
      names.push_back(name);
  }
  std::sort(names.begin(), names.end());

  if (names.empty()) {
    Error("plot_hpge_scan_overlay", "No %s*.root histograms found in %s", prefix, inputDir);
    return;
  }

  TCanvas c("c_hpge_scan_overlay", "HPGe deposited energy scan", 1100, 750);
  c.SetLogy();

  TLegend leg(0.68, 0.62, 0.90, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);

  std::vector<TFile*> openFiles;
  std::vector<TH1D*> hists;
  const int colors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kBlack};

  double maxY = 0;
  for (size_t i = 0; i < names.size(); ++i) {
    TString path = TString::Format("%s/%s", inputDir, names[i].Data());
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie())
      continue;
    TH1D* h = (TH1D*)f->Get("h_hpge_energy");
    if (!h)
      continue;
    h->SetDirectory(0);
    h->SetLineColor(colors[i % (sizeof(colors) / sizeof(colors[0]))]);
    h->SetLineWidth(2);
    hists.push_back(h);
    openFiles.push_back(f);
    maxY = std::max(maxY, h->GetMaximum());

    TString label = names[i];
    label.ReplaceAll(prefix, "Co-60 x=");
    label.ReplaceAll(".root", " mm");
    leg.AddEntry(h, label, "l");
  }

  if (hists.empty()) {
    Error("plot_hpge_scan_overlay", "No h_hpge_energy histograms found");
    return;
  }

  hists[0]->SetTitle("HPGe deposited energy scan;E_{dep} (MeV);Counts / 2 keV");
  hists[0]->SetMaximum(std::max(1.0, maxY * 2.0));
  hists[0]->Draw("hist");
  for (size_t i = 1; i < hists.size(); ++i)
    hists[i]->Draw("hist same");
  leg.Draw();

  TString pngName = TString::Format("%s.png", outputBase);
  c.SaveAs(pngName);
}
