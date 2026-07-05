void plot_hpge_energy(const char* input = "root/sim/background_plus_co60_x-800.root",
                      const char* outputBase = "root/ana/background_plus_co60_x-800") {
  gROOT->SetBatch(kTRUE);
  gSystem->AddDynamicPath("../../../NPLib/lib");
  gSystem->Load("libNPCoaxial_Germanium");

  TFile in(input, "READ");
  if (in.IsZombie()) {
    Error("plot_hpge_energy", "Cannot open %s", input);
    return;
  }

  TTree* tree = (TTree*)in.Get("SimulatedTree");
  if (!tree) {
    Error("plot_hpge_energy", "SimulatedTree not found in %s", input);
    return;
  }

  TH1D* h = new TH1D("h_hpge_energy", "HPGe deposited energy;E_{dep} (MeV);Counts / 2 keV", 1500, 0.0, 3.0);
  const Long64_t selected = tree->Draw("Coaxial_Germanium.fCoaxial_Germanium_Energy>>h_hpge_energy", "", "goff");
  Info("plot_hpge_energy", "Filled %lld HPGe entries from %s", selected, input);

  TString rootName = TString::Format("%s.root", outputBase);
  TString pngName = TString::Format("%s.png", outputBase);

  TFile out(rootName, "RECREATE");
  h->Write();
  out.Close();

  TCanvas c("c_hpge_energy", "HPGe deposited energy", 1000, 700);
  c.SetLogy();
  h->SetLineColor(kBlue + 2);
  h->SetLineWidth(2);
  h->Draw("hist");
  c.SaveAs(pngName);
}
