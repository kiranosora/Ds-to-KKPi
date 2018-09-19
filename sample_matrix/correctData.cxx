#include "utility.h"


void correctData(int recal=0){
    /**
     * 0. Get Bkg hist from generic MC
     * */
    TFile *fM =new TFile(fSig);
    TTree* tBkg=(TTree*)fM->Get("DsDecay");
    TH1F *m12Bkg=new TH1F("m12Bkg","m12Bkg",100, mBeg, mEnd);
    char cut[500];
    sprintf(cut,"costheta>-1&&has10221==0&&isSig==0&&m12 > %f && m12 <%f",mBeg, mEnd);
    tBkg->Project(m12Bkg->GetName(), "m12",cut);
    
    /**
     * 1. Get data hist after corrected
     * */
    fR=new TFile("../sample_matrix_result/corrected.root","recreate");
    binWidth = (mEnd -mBeg)/Nbin;
    binM13 = (m13End -m13Beg)/Nbin;
    binCos = (cosEnd -cosBeg)/Nbin;
    TFile *fD =new TFile(fDName);
    //TFile *fD =new TFile("../../data_703/alg_revised/ana/result.root");
    TTree* tD0 =(TTree*) fD->Get("truth");

    if(!tD0){
        cout<<" No 'truth' "<<endl;
        tD0=(TTree*) fD->Get("DsDecay"); 
    }
    else{
        isGenerator=true;
    }
    if(!tD0){
        cout<<"ERROR: NULL tree for data"<<endl;
        _exit(0);
    }

    //sprintf(cut,"m12 > %f && m12 <%f",mBeg, mEnd);
    //sprintf(cut,"isSig==1&&m12 > %f && m12 <%f",mBeg, mEnd);
    //sprintf(cut,"match==1&&m12 > %f && m12 <%f",mBeg, mEnd);
    sprintf(cut,"isSig==1&&tagMatch==1&&m12 > %f && m12 <%f",mBeg, mEnd);
    cout<<"cut: "<<cut<<endl;
    fR->cd();
    TTree* tD=tD0->CopyTree(cut);
    bool isBkg(false);
    double costheta(0), m12(0), m13(0), s12(0);
    tD->SetBranchAddress("hel1213",&costheta);
    tD->SetBranchAddress("m12",     &m12);
    tD->SetBranchAddress("s12",     &s12);
    tD->SetBranchAddress("m13",     &m13);
    
    TH1F *m12H=new TH1F("m12H","m12H",Nbin, mBeg, mEnd);
    tD->Project(m12H->GetName(), "m12",cut);
    int Nbkg=int(m12H->GetEntries() * bkgRatio);
    cout<<"  cut:  "<<cut<<endl;
    cout<<" Nbkg = "<<Nbkg<<endl;



    TFile *fMD;
    double *eff0Array=new double[Nbin*Nbin];
    double *effArray=new double[Nbin*Nbin];
    double *eff0Array_Inverse=new double[Nbin*Nbin];
    double *Ni =new double[Nbin*Nbin];
    double *NGen =new double[Nbin*Nbin];
    double *xM13 =new double[Nbin*Nbin];
    double *yM12 =new double[Nbin*Nbin];
    ifstream iEff("../sample_matrix_result/effArray.dat");
    //string astr;
    //iEff >> astr >> astr >> astr >> astr >> astr;
    for(int iN=0; iN < Nbin*Nbin; iN++){
        //iEff >> xM13[iN] >> yM12[iN] >> Ni[iN] >> NGen[iN] >> eff0Array[iN];
       iEff >>  effArray[iN];
       cout<<" effArray["<<iN<<"]="<< effArray[iN]<<endl;
    }
    iEff.close();
    //double *effArray=new double[Nbin];
    double *effArray_Inverse=new double[Nbin];
    double *xm12 = new double[Nbin];
    double *x12 = new double[Nbin];
    double *cos12 = new double[Nbin];
    TH1F *m12eff_0=new TH1F("m12eff_0","m12eff_0",Nbin, mBeg, mEnd);
    for(int binIdx=0; binIdx < Nbin; binIdx++){
        xm12[binIdx] = mBeg + (binIdx+0.5) * binWidth;
        x12[binIdx] = xm12[binIdx]* xm12[binIdx];
        cos12[binIdx] = cosBeg + (binIdx+0.5)*binCos;
        std::cout<<"x12["<<binIdx<<"] = "<<x12[binIdx]<<std::endl;
    }
    fR->cd();
    gStyle->SetMarkerStyle(2.0);
    gStyle->SetTitleSize(0.045);
    TGraph2D *gr0 =new TGraph2D(Nbin*Nbin, xM13, yM12, eff0Array);
    gr0->SetMarkerColor(kRed);
    gr0->Write("eff0");
    TGraph2D *grNSel =new TGraph2D(Nbin*Nbin, xM13, yM12, Ni);
    grNSel->SetMarkerColor(kRed);
    grNSel->Write("NSel");
    TGraph2D *grNGen =new TGraph2D(Nbin*Nbin, xM13, yM12,  NGen);
    cout<<"NGen[0]="<<NGen[0]<<endl;
    grNGen->SetMarkerColor(kRed);
    grNGen->Write("NGen");
    
    /**
     * 2. Get yk, k=0,1,2
     * */
    gSystem->Load("libMathMore");
    TF1 *L0=new TF1("L0","ROOT::Math::legendre([0],x)",-1,1);
    L0->SetParameters(0,0.0);
    TF1 *L1=new TF1("L1","ROOT::Math::legendre([0],x)",-1,1);
    L1->SetParameters(1,0.0);
    TF1 *L2=new TF1("L2","ROOT::Math::legendre([0],x)",-1,1);
    L2->SetParameters(2,0.0);

    int binIdx=0;
    double *y0=new double[Nbin];
    double *y1=new double[Nbin];
    double *y2=new double[Nbin];
    double *sWave =new double[Nbin]; 
    double *pWave =new double[Nbin]; 
    double *y0_bkg=new double[Nbin];
    double *y1_bkg=new double[Nbin];
    double *y2_bkg=new double[Nbin];
    double *sWave_bkg =new double[Nbin]; 
    double *pWave_bkg =new double[Nbin]; 
    double *delPhiP =new double[Nbin]; 
    double *delPhiM =new double[Nbin]; 
    double *dY=new double[Nbin];
    TH1F* hSWave = new TH1F("hSWave","|S|^{2}",Nbin, hBeg, hEnd);
    TH1F* hPWave = new TH1F("hPWave","|P|^{2}",Nbin, hBeg, hEnd);
    TH1F* hY0 = new TH1F("hY0", "<Y_{0}^{0}>",Nbin, hBeg, hEnd);
    TH1F* hY1 = new TH1F("hY1", "<Y_{1}^{0}>",Nbin, hBeg, hEnd);
    TH1F* hY2 = new TH1F("hY2", "<Y_{2}^{0}>",Nbin, hBeg, hEnd);
    TH1F* hY0_bkg = new TH1F("hY0_bkg", "<Y_{0}^{0}>_bkg",Nbin, hBeg, hEnd);
    TH1F* hY1_bkg = new TH1F("hY1_bkg", "<Y_{1}^{0}>_bkg",Nbin, hBeg, hEnd);
    TH1F* hY2_bkg = new TH1F("hY2_bkg", "<Y_{2}^{0}>_bkg",Nbin, hBeg, hEnd);
    TH1F* cosHel = new TH1F("cosHel", "cosHel",Nbin, -1, 1);
    TH1F* cosHel0 = new TH1F("cosHel0", "cosHel0",Nbin, -1, 1);
    for(int binIdx=0; binIdx < Nbin; binIdx++){
        y0[binIdx]=0;
        y1[binIdx]=0;
        y2[binIdx]=0;
        sWave[binIdx]=0;
        pWave[binIdx]=0;
        y0_bkg[binIdx]=0;
        y1_bkg[binIdx]=0;
        y2_bkg[binIdx]=0;
        sWave_bkg[binIdx]=0;
        pWave_bkg[binIdx]=0;
        delPhiP[binIdx]=0;
        delPhiM[binIdx]=0;
        dY[binIdx] =0;
    }

    cout<<setiosflags(ios::floatfield);
    double pi = TMath::Pi();
    /**
     * a. cal intergration of data
     * */

    double phspFactor=0.;
    int cosIdx=0;
    double eff=0;
    for(int iD = 0; iD < tD->GetEntries(); iD++){
        tD->GetEntry(iD);
        if(fabs(costheta)>1) {
            cout<<"iD= "<<iD<<" costheta="<<costheta<<endl;
            continue;
        }
        binIdx =int ((m12 - mBeg)/binWidth);
        cosIdx =int ((costheta - cosBeg)/binCos);
        if(binIdx < 0 || binIdx >= Nbin) continue;
        if(cosIdx < 0 || cosIdx >= Nbin) continue;
        if(binIdx==Nbin){
            binIdx--;
        }

        if(binIdx < 0 || binIdx >= Nbin){
            cout<<"m12= "<<m12<<" mBeg= "<<mBeg<<" del= "<< m12-mBeg<<" binWidth= "<<binWidth<<" binIdx = " << binIdx<<endl;
        }
        eff = effArray[cosIdx * Nbin + binIdx];
        //eff = 1/eff0Array[cosIdx * Nbin + binIdx];
        //eff=1;
        if(isGenerator) eff=1;
        double t=sqrt(5/(4.0*pi))* L2->Eval(costheta) * eff;
        phspFactor = sqrt( 1 - 4*mK*mK / s12);
        //phspFactor = 1.0;
        if(eff<0){
            cout<<" error: pol4->Eval("<<costheta<<")= "<< eff <<endl;
        }
        m12eff_0->Fill(m12, eff);
        //cout<<" m12 = "<<m12<<" eff = "<<eff<<" effArray["<<cosIdx<<"*"<<Nbin<<"+"<<binIdx<<"]"<<endl;
        //return;
        y0[binIdx] += sqrt(1/(4.0*pi))* L0->Eval(costheta) * eff / phspFactor;
        hY0->Fill(m12,sqrt(1/(4.0*pi))* L0->Eval(costheta) * eff/ phspFactor); 
        y1[binIdx] += sqrt(3/(4.0*pi))* L1->Eval(costheta) * eff/ phspFactor;
        hY1->Fill(m12,sqrt(3/(4.0*pi))* L1->Eval(costheta) * eff/ phspFactor); 
        y2[binIdx] += sqrt(5/(4.0*pi))* L2->Eval(costheta) * eff/ phspFactor;
        hY2->Fill(m12,sqrt(5/(4.0*pi))* L2->Eval(costheta) * eff/ phspFactor); 
        cosHel->Fill(costheta, eff); 
        cosHel0->Fill(costheta); 
        dY[binIdx] += sqrt(4*pi)*y0[binIdx] - sqrt(5*pi)*y2[binIdx];
        //cout<<" y2["<<binIdx<<" ]= "<<y2[binIdx]<<" L2->Eval(costheta)= "<<L2->Eval(costheta)<<" UseEff["<<binIdx<<"] = "<<eff<<endl;
    }


    /**
     * b. cal intergration of bkg( from generic MC)
     * */
    double norm = Nbkg*1.0 / tBkg->GetEntries();
    norm=0;
    TH1F *m12eff_bkg=new TH1F("m12eff_bkg","m12eff_bkg",Nbin, mBeg, mEnd);
    tBkg->SetBranchAddress("hel1213",&costheta);
    tBkg->SetBranchAddress("m12",     &m12);
    tBkg->SetBranchAddress("s12",     &s12);
    tBkg->SetBranchAddress("m13",     &m13);
    for(int iBkg = 0; iBkg < tBkg->GetEntries(); iBkg++){
        tBkg->GetEntry(iBkg);
        binIdx =int ((m12 - mBeg)/binWidth);
        cosIdx =int ((costheta - cosBeg)/binCos);
        if(fabs(costheta)>1) {
            cout<<"costheta="<<costheta<<endl;
            _exit(0);
        }
        if(binIdx==Nbin){
            binIdx--;
        }
        if(binIdx < 0 || binIdx >= Nbin){
            cout<<"m12= "<<m12<<" mBeg= "<<mBeg<<" del= "<< m12-mBeg<<" binWidth= "<<binWidth<<" binIdx = " << binIdx<<endl;
        }
        eff = effArray[cosIdx * Nbin + binIdx];
        //eff = 1/eff0Array[cosIdx * Nbin + binIdx];
        //eff = 1;
        if(isGenerator) eff=1;
        double t=sqrt(5/(4.0*pi))* L2->Eval(costheta) * eff;
        if(eff <0){
            cout<<" error: UseEff["<<binIdx<<"]= "<<eff<<endl;
        }
        phspFactor = sqrt( 1 - 4*mK*mK / s12);
        //phspFactor = 1.0;
        m12eff_bkg->Fill(m12, norm*eff/ phspFactor);
        y0_bkg[binIdx] += sqrt(1/(4.0*pi))* L0->Eval(costheta) * eff/ phspFactor;
        hY0_bkg->Fill(m12,sqrt(1/(4.0*pi))* L0->Eval(costheta) * eff/ phspFactor); 
        y1_bkg[binIdx] += sqrt(3/(4.0*pi))* L1->Eval(costheta) * eff/ phspFactor;
        hY1_bkg->Fill(m12,sqrt(3/(4.0*pi))* L1->Eval(costheta) * eff/ phspFactor); 
        y2_bkg[binIdx] += sqrt(5/(4.0*pi))* L2->Eval(costheta) * eff/ phspFactor;
        hY2_bkg->Fill(m12,sqrt(5/(4.0*pi))* L2->Eval(costheta) * eff/ phspFactor); 
    }

    /**
     * c. cal intergration of signal(data)
     * */
    cout<<" norm = "<<Nbkg<<" / "<<tBkg->GetEntries()<<" = "<<norm<<endl;
    //norm=0;
    for(int bI=0; bI<Nbin; bI++){
        y0[bI]  = y0[bI] - norm*y0_bkg[bI];
        y1[bI]  = y1[bI] - norm*y1_bkg[bI];
        double a = y2[bI];
        y2[bI]  = y2[bI] - norm*y2_bkg[bI];

        cout<<" y0["<<bI<<"] = "<<y0[bI]<<" y2["<<bI<<"] = "<<a<<" - "<< norm<<" * "<< y2_bkg[bI] << " = "<< y2[bI]<<" y0["<<bI<<"]/y2["<<bI<<"="<<y0[bI]/y2[bI]<<endl;
    }

    /**
     * 3. cal S-Wave, P-Wave and phase difference
     * 
     * |P|= sqrt (sqrt(5.0)/2.0 * sqrt(4*pi) * y2);
     * |S| = sqrt( sqrt(4*pi) * y0 - |P|*|P|);
     * cos(delPhi) = sqrt(4*TMath::pi) *y1 /( 2* |S| *|P|);
     * */
    double t1=0;
    double nS=0;
    double nP=0;
    for(int binIdx=0; binIdx<Nbin; binIdx++){
        t1=sqrt(5.0*pi) * y2[binIdx];
        //t1=sqrt(5.0)/2.0 * sqrt(4*pi) * y2[binIdx];
        double t2=sqrt(4*pi) * y0[binIdx]  - sqrt(5*pi)*y2[binIdx];
        hSWave->SetBinContent(binIdx+1,sqrt(4*pi) * y0[binIdx]  - sqrt(5*pi)*y2[binIdx]);
        hPWave->SetBinContent(binIdx+1, sqrt(5*pi)*y2[binIdx]);
        nS += sqrt(4*pi) * y0[binIdx]  - sqrt(5*pi)*y2[binIdx];
        nP += sqrt(5*pi)*y2[binIdx];
        if(y2[binIdx]<0){
            cout<<"y2["<<binIdx<<"] = "<<y2[binIdx]<<endl;
            y2[binIdx]=0;
            t2=0;
            //t1=0;
            //exit(0);
        }
        pWave[binIdx] = sqrt(t1);   
        sWave[binIdx] = sqrt(t2 );   
        delPhiP[binIdx] = acos( sqrt(4*pi) * y1[binIdx]/(2* sWave[binIdx] * pWave[binIdx])); 
        delPhiM[binIdx] = -acos( sqrt(4*pi) * y1[binIdx]/(2* sWave[binIdx] * pWave[binIdx])); 
    }
    cout<<" nS / nP = "<< nS<<" / "<< nP<< " = "<< nS/nP<<endl;
    //return;
    TH1F *m12eff=new TH1F((*m12eff_0) - (*m12eff_bkg));
    m12eff->Write();
    m12H->Write();
    cosHel->Write();
    cosHel0->Write();
    hSWave->Draw();
    hSWave->GetXaxis()->SetNdivisions(505);
    hSWave->GetXaxis()->SetTitle("m(K^{+}K^{-})(GeV/c^{2})");
    hSWave->GetXaxis()->CenterTitle();
    hSWave->Write(); 
    hPWave->Draw();
    hPWave->GetXaxis()->SetNdivisions(505);
    hPWave->GetXaxis()->SetTitle("m(K^{+}K^{-})(GeV/c^{2})");
    hPWave->GetXaxis()->CenterTitle();
    hPWave->Write(); 
    TGraph * gP = new TGraph(Nbin, xm12, pWave);
    gP->Write("pWave");
    TGraph * gS = new TGraph(Nbin, xm12, sWave);
    gS->Write("sWave");
    TGraph * gPhiP = new TGraph(Nbin, xm12, delPhiP);
    gPhiP->Write("delPhiP");
    TGraph * gPhiM = new TGraph(Nbin, xm12, delPhiM);
    gPhiM->Write("delPhiM");
    TGraph *gY0 = new TGraph(Nbin, xm12, y0);
    gY0->Draw("ap");
    gY0->Write("Y0");
    TGraph *gY1 = new TGraph(Nbin, xm12, y1);
    gY1->Draw("ap");
    gY1->Write("Y1");
    TGraph *gY2 = new TGraph(Nbin, xm12, y2);
    gY2->Draw("ap");
    gY2->Write("Y2");
    TGraph *gdY = new TGraph(Nbin, xm12, dY);
    gdY->Draw("ap");
    gdY->Write("dY");
    hY0_bkg->SetMarkerStyle(20);
    //hY0->SetMarkerSize(2);
    hY0_bkg->Draw("E");
    hY0_bkg->GetXaxis()->SetNdivisions(505);
    hY0_bkg->GetXaxis()->SetTitle("m(K^{+}K^{-})(GeV/c^{2})");
    hY0_bkg->GetXaxis()->CenterTitle();
    hY0_bkg->Write();
    hY1_bkg->SetMarkerStyle(20);
    //hY1->SetMarkerSize(2);
    hY1_bkg->Draw("E");
    hY1_bkg->GetXaxis()->SetNdivisions(505);
    hY1_bkg->GetXaxis()->SetTitle("m(K^{+}K^{-})(GeV/c^{2})");
    hY1_bkg->GetXaxis()->CenterTitle();
    hY1_bkg->Write();
    hY2_bkg->SetMarkerStyle(20);
    //hY2->SetMarkerSize(2);
    hY2_bkg->Draw("E");
    hY2_bkg->GetXaxis()->SetNdivisions(505);
    hY2_bkg->GetXaxis()->SetTitle("m(K^{+}K^{-})(GeV/c^{2})");
    hY2_bkg->GetXaxis()->CenterTitle();
    hY2_bkg->Write();
    TCanvas *cc=new TCanvas("cc","cc",1600,600);
    cc->Divide(3,1);
    cc->cd(1);
    hY0->SetMarkerStyle(20);
    //hY0->SetMarkerSize(2);
    hY0->Draw("E");
    hY0->GetXaxis()->SetNdivisions(505);
    hY0->GetXaxis()->SetTitle("m(K^{+}K^{-})(GeV/c^{2})");
    hY0->GetXaxis()->CenterTitle();
    hY0->Write();
    cc->cd(2);
    hY1->SetMarkerStyle(20);
    //hY1->SetMarkerSize(2);
    hY1->Draw("E");
    hY1->GetXaxis()->SetNdivisions(505);
    hY1->GetXaxis()->SetTitle("m(K^{+}K^{-})(GeV/c^{2})");
    //hY1->GetXaxis()->SetTitleCenter();
    hY1->GetXaxis()->CenterTitle();
    hY1->Write();
    cc->cd(3);
    hY2->SetMarkerStyle(20);
    //hY2->SetMarkerSize(2);
    hY2->Draw("E");
    hY2->GetXaxis()->SetNdivisions(505);
    hY2->GetXaxis()->SetTitle("m(K^{+}K^{-})(GeV/c^{2})");
    //hY2->GetXaxis()->SetTitleCenter();
    hY2->GetXaxis()->CenterTitle();
    hY2->Write();
    delete[] NGen;
}


int main(int argc, char **argv){
    //fR=new TFile("corrected.root","recreate");
    //binWidth = (mEnd -mBeg)/Nbin;
    //calEff();
    if(argc > 1){
        correctData(atoi(argv[1]));
    }
    else{
        correctData();
    }
}
