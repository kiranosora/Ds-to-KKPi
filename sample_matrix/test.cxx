#include "utility.h"

double * calEff(double* eff0Array, double **Ni, double *NGen, double *xM13, double *yM12, TMatrixD *inverseMEff){
    char cut[500];
    double genM13Beg=0;
    double genM13End=0;
    double genM12Beg=0;
    double genM12End=0;
    double selM13Beg=0;
    double selM13End=0;
    double selM12Beg=0;
    double selM12End=0;

    TChain *tc =new TChain("truth","truth");
    tc->AddFile(fGen);
    for(int iGen=0; iGen < Nbin; iGen++){
        genM13Beg = cosBeg + binCos * iGen;
        genM13End = cosBeg + binCos * (iGen+1);
        for(int iGenM=0; iGenM< Nbin; iGenM++){
            genM12Beg = mBeg + binWidth * iGenM;
            genM12End = mBeg + binWidth * (iGenM+1);
            xM13[iGen*Nbin+iGenM] = genM13Beg + 0.5 * binM13;
            yM12[iGen*Nbin+iGenM] = genM12Beg + 0.5 * binWidth;
            cout<<" xM13["<<iGen*Nbin+iGenM<<"]="<<xM13[iGen*Nbin+iGenM]<<" genM13Beg["<<iGen*Nbin+iGenM<<"]="<<genM13Beg<<endl;
        }
    }


    /**
     * 1. Get Entries at generator level
     * */
    double  m12(0), m13(0), costheta(0);
    tc->SetBranchAddress("hel1213",     &costheta);
    tc->SetBranchAddress("m12",     &m12);
    int m13Idx(0), m12Idx(0);
    cout<<" read generator level truth "<<endl;
    for(int itc=0; itc< tc->GetEntries(); itc++){
        tc->GetEntry(itc);
        if(costheta < cosBeg || costheta > cosEnd) continue;
        if(m12 < mBeg || m12 > mEnd) continue;
        m13Idx =( costheta - cosBeg)/binCos;
        m12Idx =( m12 - mBeg)/binWidth;
        NGen[m13Idx*Nbin+m12Idx]++;
    }
    
    for(int iN=0; iN< Nbin*Nbin; iN++){
        cout<<__LINE__<<" NGen["<<iN<<"]="<<NGen[iN]<<endl;
    }

    /**
     * 2. Get Entries after selection and efficiency matrix
     * */
    TH1F *hGen=new TH1F("hGen","hGen",100, mBeg, mEnd);
    //TFile *f =new TFile("../alg_revised/ChooseBestSigTagPair/root/0001_abs.root");
    TFile *fM =new TFile(fPhsp);
    TFile *tem=new TFile("result/tmp.root","recreate"); 
    TTree* t0=(TTree*) fM->Get("DsDecay");
    TFile *fT =new TFile("result/tmp.root","recreate");
    TTree* t=t0->CopyTree(cut);
    //TTree* t=(TTree*) fM->Get("KKPi");
    int *Nij= new int[Nbin*Nbin];
    sprintf(cut,"isSig==1&&tagMatch==1&&m12 > %f && m12 <%f",mBeg, mEnd);
    int Nent= t->Project("hGen","m12",cut);
    hGen->GetEntries();;
    cout<<" Total Entries used: "<< Nent<< endl;
    cout<<"cut: "<<cut<<endl;
    //hGen->Write("hGen");
    double m12_truth(0), costheta_truth(0);
    t->SetBranchAddress("hel1213",     &costheta);
    t->SetBranchAddress("m12",     &m12);
    t->SetBranchAddress("hel1213_truth",     &costheta_truth);
    t->SetBranchAddress("m12_truth",     &m12_truth);
    int m13IdxTruth(0), m12IdxTruth(0);
    cout<<" read info after selection "<<endl;
    for(int it=0; it< t->GetEntries(); it++){
        t->GetEntry(it);
        if(costheta < cosBeg || costheta > cosEnd) continue;
        if(m12 < mBeg || m12 > mEnd) continue;
        if(costheta_truth < cosBeg || costheta_truth > cosEnd) continue;
        if(m12_truth < mBeg || m12_truth > mEnd) continue;
        if(it%1000==0) cout<<" "<<it<<" events processed"<<endl;
        m13Idx =( costheta - cosBeg)/binCos;
        m13IdxTruth = ( costheta_truth - cosBeg)/binCos;
        m12Idx =( m12 - mBeg)/binWidth;
        m12IdxTruth = ( m12_truth - mBeg)/binWidth;
        if(m13Idx<0|| m13Idx>=Nbin || m12Idx<0 || m12Idx>=Nbin){
            cout<<" costheta = "<<costheta<<" m12 = "<<m12<<endl;
            continue;
        }
        Ni[m13Idx*Nbin+m12Idx][m13IdxTruth*Nbin + m12IdxTruth]++;
        //Nij[m13Idx*Nbin+m12Idx]++;
        if(it%1000==0) cout<<__LINE__<<" Ni["<<m13Idx*Nbin+m12Idx<<"]["<<m13IdxTruth*Nbin + m12IdxTruth<<"]="<<Ni[m13Idx*Nbin+m12Idx][m13IdxTruth*Nbin + m12IdxTruth] <<" events processed"<<endl;
    }
    double *eff = new double[Nbin * Nbin];
    for(int rowIdx = 0; rowIdx < Nbin*Nbin; rowIdx++){
        cout<<"Nij["<<rowIdx<<"]="<<Nij[rowIdx]<<endl;
        for(int rowIdxM=0; rowIdxM < Nbin*Nbin; rowIdxM++){
            cout<<__LINE__<<" Ni["<<rowIdx<<"]["<<rowIdxM<<"]= "<<Ni[rowIdx][rowIdxM]<<" NGen["<<rowIdxM<<"]="<<NGen[rowIdxM]<<endl;
            eff0Array[rowIdx*Nbin*Nbin+rowIdxM]= Ni[rowIdx][rowIdxM]*1.0/NGen[rowIdxM];
            cout<<__LINE__<<" eff0Array["<<rowIdx*Nbin*Nbin+rowIdxM<<"] = "<<eff0Array[rowIdx*Nbin*Nbin+rowIdxM]<<endl;
            //cout<<"      xM13["<<rowIdx*Nbin + rowIdxM<<"] = "<<xM13[rowIdx*Nbin + rowIdxM]<<endl;
            //cout<<"      yM12["<<rowIdx*Nbin + rowIdxM<<"] = "<<yM12[rowIdx*Nbin + rowIdxM]<<endl;
            //cout<<" eff0Array["<<rowIdx*Nbin + rowIdxM<<"] = "<<eff0Array[rowIdx*Nbin + rowIdxM]<<endl;
            if(isGenerator){
                cout<<"isGenerator: "<<isGenerator<<" Ni["<<rowIdx<<"]= "<<Ni[rowIdx]<<" cut: "<<cut<<endl;
                if( rowIdx == rowIdxM){
                    eff0Array[rowIdx*Nbin*Nbin+rowIdxM]=1;
                }
                else{
                    eff0Array[rowIdx*Nbin*Nbin+rowIdxM]=0;
                }
            }
        }
    }

    double *effArray = new double[Nbin*Nbin];

    delete tc;
    delete t;
    delete[] Nij;
    delete[] eff;
    delete tem;
    return effArray;
    //return eff0Array ;
}

void correctData(int recal=0){
    char cut[500];
    binWidth = (mEnd -mBeg)/Nbin;
    binM13 = (m13End -m13Beg)/Nbin;
    binCos = (cosEnd -cosBeg)/Nbin;


    double *eff0Array=new double[Nbin*Nbin*Nbin*Nbin];
    std::uninitialized_fill_n(eff0Array, Nbin*Nbin*Nbin*Nbin, 0);
    double **Ni  =new double*[Nbin*Nbin];
    for(int ii =0; ii < Nbin*Nbin; ii++){
        Ni[ii] = new double[Nbin*Nbin];
        std::uninitialized_fill_n(Ni[ii], Nbin*Nbin, 0);
    }
    //int *Ni =new int[Nbin*Nbin];
    //std::uninitialized_fill_n(Ni, Nbin*Nbin, 0);
    double *NGen =new double[Nbin*Nbin];
    std::uninitialized_fill_n(NGen, Nbin*Nbin, 0);
    double *xM13 =new double[Nbin*Nbin];
    double *yM12 =new double[Nbin*Nbin];
    TMatrixD *inverseMEff;
    double *effArray=calEff(eff0Array, Ni, NGen, xM13, yM12, inverseMEff);
    ofstream ofs("result/array.dat");
    ofs<<"\t xM13\t\t yM12\t\t Ni\t\t NGen\t\t eff0Array\t"<<endl;
    for(int i =0; i< Nbin*Nbin; i++){
        ofs<<"\t "<<xM13[i]<<"\t "<<yM12[i]<<"\t\t "<<Ni[i]<<"\t "<<NGen[i]<<"\t"<<effArray[i]<<endl;
    }
    ofs.close();

}


int main(int argc, char **argv){
    if(argc > 1){
        correctData(atoi(argv[1]));
    }
    else{
        correctData();
    }
}
