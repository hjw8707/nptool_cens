TTree* LoadTree(const char* type, const char* exst);
TH1* GetSpectrum(const char* cut = NULL, Double_t upp = 3, Int_t nBin = 300,
		 Bool_t flagDC = true, Bool_t flagADD = true,Int_t iEntry = 0);
TFitResult* FitSpectrum(TH1* h1, Double_t fitPos, Double_t range = 0);

const char* gcut = "ELab > 35";
TFile *file = NULL;
TTree *tree = NULL;
TString type;
TString exst;
TDirectory *tdir = new TDirectory("tdir","tdir");

TTree* LoadTree(const char* _type, const char* _exst) {
  file = new TFile(Form("Analysis/Sr78_ana_%s_%s.root",_exst,_type));
  tree = static_cast<TTree*>(file->Get("PhysicsTree"));
  type = _type;
  if      (type.BeginsWith("dali")) type = "dali";
  else if (type.BeginsWith("csi"))  type = "csi";
  else                              type = "grape";
  exst = _exst;
  return tree;}

TH1* GetSpectrum(const char* _cut, Double_t upp, Int_t nBin,
		 Bool_t flagDC, Bool_t flagADD, Int_t iEntry) {
  if (!tree) return NULL;

  TH1 *h1;
  TString brName(type.Data());
  brName += "E";
  if (flagDC) brName += "dc";
  if (flagADD) brName += "add";
  TString hName("h");
  hName += brName;
  hName += "_";
  hName += exst;
  tree->Draw(Form("%s>>%s(%d,0,%f)",brName.Data(),hName.Data(),nBin,upp),_cut);
  h1 = static_cast<TH1*>(gDirectory->Get(hName.Data()));
  h1->SetTitle(";Energy [MeV];Event #");
  h1->SetDirectory(tdir);
  return h1;}

TFitResult* FitSpectrum(TH1* h1, Double_t fitPos, Double_t range) {
  if (!h1) return NULL;

  if (range <= 0) range = 0.08;
  TFitResultPtr frp = h1->Fit("gausn","SQ",0,fitPos*(1 - range), fitPos*(1 + range));
  return frp.Get();}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void LoadAllExGammas(const char* type, const char* cut = gcut) {
  TCanvas *c1 = new TCanvas;
  c1->Divide(2,2);
  const char* exs[] = { "1st0", "1st2", "2nd0",  "2nd2" };
  Double_t pos_mcsm[4] = { 0., 0.289, 2.251, 0.701 };
  Double_t pos_bmf[4] = { 0., 0.289, 0.471, 0.811 };
  Double_t *pos;
  TString ttype(type);
  if (ttype.Contains("mcsm")) pos = pos_mcsm;
  else                        pos = pos_bmf;
  
  for (int i = 0 ; i < 4 ; i++) {
    if (file) { file->Close(); delete file; }
    c1->cd(i+1);
    LoadTree(type, exs[i]);
    TH1 *h1 = GetSpectrum(cut);
    if (i != 0) {
      TFitResult* tfr = FitSpectrum(h1, pos[i]);
      Double_t nEvFit = tfr->GetParams()[0] / h1->GetBinWidth(1);
      Double_t nTotEv = double(tree->GetEntries());
      std::cout << "Eff. for " << h1->GetName() << ": " << (nEvFit/nTotEv) << std::endl;}
    h1->Draw();}
  tdir->cd();}
