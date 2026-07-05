void plot_gas_energy(const char* input = "root/sim/cssu_cf4_alpha.root",
                     const char* outputBase = "root/ana/cssu_cf4_alpha") {
  gROOT->SetBatch(kTRUE);
  gSystem->AddDynamicPath("../../../NPLib/lib");
  gSystem->Load("libNPCSSU");

  TFile in(input, "READ");
  if (in.IsZombie()) {
    Error("plot_gas_energy", "Cannot open %s", input);
    return;
  }

  TTree* tree = (TTree*)in.Get("SimulatedTree");
  if (!tree) {
    Error("plot_gas_energy", "SimulatedTree not found in %s", input);
    return;
  }

  TH1D* h = new TH1D("h_cssu_gas_energy", "CSSU CF4 gas deposited energy;E_{dep} (MeV);Counts / 5 keV",
                     1200, 0.0, 6.0);
  Long64_t selected = tree->Draw("CSSU.fCSSU_Energy>>h_cssu_gas_energy", "", "goff");
  Info("plot_gas_energy", "Filled %lld CSSU gas entries from %s", selected, input);

  TString rootName = TString::Format("%s.root", outputBase);
  TString pngName = TString::Format("%s.png", outputBase);
  TFile out(rootName, "RECREATE");
  h->Write();
  out.Close();

  TCanvas c("c_cssu_gas_energy", "CSSU CF4 gas deposited energy", 1000, 700);
  h->SetLineColor(kBlue + 2);
  h->SetLineWidth(2);
  h->Draw("hist");
  c.SaveAs(pngName);
}
