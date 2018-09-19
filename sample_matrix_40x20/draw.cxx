void draw(){
    gStyle->SetOptStat(0);
    TFile *f=new TFile("../sample_matrix_result/corrected.root");
    TLine *line=new TLine(0.95, 0, 1.15,0);
    line->SetLineStyle(2);
    
    
    TCanvas *cc=new TCanvas("cc","cc",899*3*2.0/3,799*2.0/3);
    cc->Divide(3,1);
    
    //TCanvas *c0=new TCanvas("c0","c0",899,799);
    cc->cd(1);
    TH1F* hY0=(TH1F*)f->Get("hY0");
    hY0->Draw("E");
    line->Draw();
    hY0->GetXaxis()->SetTitleSize(0.045);

    //TCanvas *c1=new TCanvas("c1","c1",899,799);
    cc->cd(2);
    TH1F* hY1=(TH1F*)f->Get("hY1");
    hY1->Draw("E");
    line->Draw();
    hY1->GetXaxis()->SetTitleSize(0.045);
    
    //TCanvas *c2=new TCanvas("c2","c2",899,799);
    cc->cd(3);
    TH1F* hY2=(TH1F*)f->Get("hY2");
    hY2->Draw("E");
    line->Draw();
    hY2->GetXaxis()->SetTitleSize(0.045);
    cc->Print("../sample_matrix_result/IndependentModel.pdf");
    cc->Print("../sample_matrix_result/Y0.png");

    TCanvas *ch=new TCanvas("ch","ch",899*2,799);
    ch->Divide(2,1);
    
    //TCanvas *c0=new TCanvas("c0","c0",899,799);
    ch->cd(1);
    TH1F* hSWave=(TH1F*)f->Get("hSWave");
    hSWave->Draw();
    cout<<" S Wave: "<<hSWave->GetSumOfWeights()<<" Nbin="<<hSWave->GetNbinsX()<<endl;
    line->Draw();
    hSWave->GetXaxis()->SetTitleSize(0.045);

    ch->cd(2);
    TH1F* hPWave=(TH1F*)f->Get("hPWave");
    hPWave->Draw();
    line->Draw();
    cout<<" P Wave: "<<hPWave->GetSumOfWeights()<<endl;
    hPWave->GetXaxis()->SetTitleSize(0.045);
    int Nbin=40;
    double nS=0;
    double nP=0;
    for(int i=0; i <= Nbin; i++){
       nS += hSWave->GetBinContent(i); 
       nP += hPWave->GetBinContent(i); 
    }
    cout<<" nS/nP = "<<nS<<"/"<<nP<<" = "<< nS*1.0/nP<<endl;
    ch->Print("../sample_matrix_result/SP.png");

    return;
    TCanvas *cEff=new TCanvas("cEff","efficiency", 900, 800);
    cEff->Divide(2,2);
    
    cEff->cd(1);
    TGraph* NGen=(TGraph*) f->Get("NGen");
    NGen->Draw("apl");
    NGen->GetXaxis()->SetTitle("m^{2}(K^{+}K^{-})");
    NGen->GetYaxis()->SetTitle("GenPHSP");
    NGen->GetXaxis()->SetTitleSize(0.045);
    NGen->GetYaxis()->SetTitleSize(0.045);
    NGen->SetTitle("GenPHSP vs m^{2}(K^{+}K^{-})");

    cEff->cd(2);
    TGraph* NSel=(TGraph*) f->Get("NSel");
    NSel->Draw("apl");
    NSel->GetXaxis()->SetTitle("m^{2}(K^{+}K^{-})");
    NSel->GetYaxis()->SetTitle("SelPHSP");
    NSel->GetXaxis()->SetTitleSize(0.045);
    NSel->GetYaxis()->SetTitleSize(0.045);
    NSel->SetTitle("SelPHSP vs m^{2}(K^{+}K^{-})");
    
    cEff->cd(3);
    TGraph* eff0=(TGraph*) f->Get("eff0");
    eff0->Draw("apl");
    eff0->GetXaxis()->SetTitle("m^{2}(K^{+}K^{-})");
    eff0->GetYaxis()->SetTitle("Efficiency");
    eff0->GetXaxis()->SetTitleSize(0.045);
    eff0->GetYaxis()->SetTitleSize(0.045);
    eff0->SetTitle("Efficiency vs m^{2}(K^{+}K^{-})");
    
    cEff->cd(4);
    TH1F* m12eff_0=(TH1F*) f->Get("m12eff_0");
    m12eff_0->Draw();
    m12eff_0->GetXaxis()->SetTitle("m(K^{+}K^{-})");
    //m12eff_0->GetYaxis()->SetTitle("Entries");
    m12eff_0->GetXaxis()->SetTitleSize(0.045);
    m12eff_0->GetYaxis()->SetTitleSize(0.045);
    m12eff_0->SetTitle("m(K^{+}K^{-}) after correction");
    cEff->Print("../sample_matrix_result/eff.png");
}
