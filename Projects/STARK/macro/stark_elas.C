{
  TFile *file = new TFile("../../Outputs/Analysis/stark_elas.root");
  TTree *tree = static_cast<TTree*>(file->Get("PhysicsTree"));

  TCanvas *c1 = new TCanvas("stark","stark",500,800);
  c1->Divide(1,2);

  TVirtualPad *vpad = c1->cd(1);
  vpad->SetLogz(false);

  tree->Draw("hPos.fZ:hPos.fY:hPos.fX>>hhitPos","","lego"); // reconstructed hit position

  vpad = c1->cd(2);
  vpad->SetLogz(false);

  tree->Draw("E:dE>>hdee(200,0,30,200,0,50)","nhit > 1","col"); // dE-E plot

}
