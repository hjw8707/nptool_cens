void draw_summary()
{
    TString data_path = "data";

    auto graph = new TGraph();
    TSystemFile *sysFile = nullptr;
    TIter nextFile(TSystemDirectory("search",data_path.Data()).GetListOfFiles());
    while ((sysFile=(TSystemFile*)nextFile()))
    {
        TString fileName = sysFile -> GetName();
        if (fileName.Index("geometrical_efficiency_")==0)
        {
            auto i1 = fileName.Index("_z");
            auto i2 = fileName.Index(".ana");
            TString z_value = fileName(i1+2,i2-i1-4);
            ifstream file(Form("%s/%s",data_path.Data(),fileName.Data()));
            double eff_value;
            cout << fileName << " = " << z_value << " " << eff_value << endl;
            file >> eff_value;
            graph -> SetPoint(graph->GetN(),z_value.Atof(),eff_value);
        }
    }
    graph -> Sort();

    auto cvs = new TCanvas("cvs","",1600,1200);
    //auto hist = new TH2D("hist",";z (mm);geometrical efficiency",100,-160,160,100,0.18,0.55);
    auto hist = new TH2D("hist",";z (mm);geometrical efficiency",100,-160,160,100,0.18,1.00);
    hist -> SetStats(0);
    hist -> Draw();
    graph -> Draw("same*l");

    gPad -> SetGridx();
    gPad -> SetGridy();
    gPad -> SaveAs("figures/efficiency.png");
    gPad -> SaveAs("figures/efficiency.pdf");
}
