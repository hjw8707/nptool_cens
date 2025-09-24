void draw_hit(TString fileName) {
    gStyle -> SetOptStat(11111111);
    cout << fileName << endl;

    TString anaName = fileName;
    anaName = anaName(anaName.Index("/")+1,anaName.Sizeof());
    anaName = anaName(0,anaName.Index(".ana.root"));
    anaName.ReplaceAll(".","_");

    auto tree = new TChain("PhysicsTree");
    tree -> AddFile(fileName);

    auto cvs = new TCanvas("cvs","",1500,600);
    cvs -> Divide(2,1);

    auto draw_ft = [cvs,tree](int i, TString expression, TCut cut="", TString option="")
    {
        TString htag = expression;
        if (htag.Index(":")>0||htag.Index("(")>0||htag.Index("[")>0) htag = Form("%d",i);
        TString expressFull = expression.Data();
        if (expressFull.Index(">>")<0) expressFull = Form("%s>>hist_%s",expression.Data(),htag.Data());
        cvs -> cd(i);
        auto entries = tree -> Draw(expressFull,cut,option);
        cout << expressFull << " " << cut.GetTitle() << " " << option << " # " << entries << endl;
    };

    TCut cut = "";

    int icvs = 1;
    draw_ft(icvs++,"nhit"  , cut);
    draw_ft(icvs++,"detN"  , cut);

    auto hist1 = (TH1D*) (gDirectory->FindObject("hist_nhit"));
    auto hist2 = (TH1D*) (gDirectory->FindObject("hist_detN"));
    double total = hist1 -> GetEntries();
    double collected = hist2 -> GetEntries();
    double efficiency = collected/total;

    TString data_path = "data";
    TString fileName2 = fileName;
    fileName2.ReplaceAll(Form("%s/",data_path.Data()),Form("%s/geometrical_efficiency_",data_path.Data()));
    fileName2.ReplaceAll(".root",".txt");
    ofstream file2(fileName2);
    cout << fileName2 << endl;
    file2 << efficiency << endl;

    cvs -> SaveAs(Form("figures/summary_%s.png",anaName.Data()));
}
