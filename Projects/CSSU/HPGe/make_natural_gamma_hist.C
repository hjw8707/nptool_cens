void make_natural_gamma_hist() {
  TFile* file = new TFile("natural_gamma_lines.root", "RECREATE");
  TH1D* hist = new TH1D("h_gamma", "Natural background gamma lines;E_{#gamma} (MeV);relative intensity", 3000, 0, 3.0);

  auto add_line = [&](double energy_keV, double intensity) {
    int bin = hist->FindBin(energy_keV / 1000.0);
    hist->SetBinContent(bin, hist->GetBinContent(bin) + intensity);
  };

  // Representative natural-background lines from K-40 and U/Th chains.
  add_line(238.632, 43.6);    // Pb-212
  add_line(295.224, 18.4);    // Pb-214
  add_line(338.320, 11.3);    // Ac-228
  add_line(351.932, 35.6);    // Pb-214
  add_line(511.000, 5.0);     // annihilation
  add_line(583.187, 85.0);    // Tl-208
  add_line(609.312, 45.5);    // Bi-214
  add_line(727.330, 6.7);     // Bi-212
  add_line(911.204, 25.8);    // Ac-228
  add_line(968.971, 15.8);    // Ac-228
  add_line(1120.287, 14.9);   // Bi-214
  add_line(1238.110, 5.8);    // Bi-214
  add_line(1377.669, 4.0);    // Bi-214
  add_line(1460.822, 10.6);   // K-40
  add_line(1764.494, 15.3);   // Bi-214
  add_line(2204.210, 5.1);    // Bi-214
  add_line(2614.511, 99.0);   // Tl-208

  hist->Write();
  file->Close();
}

