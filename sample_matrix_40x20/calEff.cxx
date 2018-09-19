#include "utility.h"
bool recal=true;
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

    int m13Idx(0), m12Idx(0);
    double  m12(0), m13(0), costheta(0);
    if(recal){
        TChain *tc =new TChain("truth","truth");
        tc->AddFile(fGen);
        for(int iGen=0; iGen < NbinCos; iGen++){
            genM13Beg = cosBeg + binCos * iGen;
            genM13End = cosBeg + binCos * (iGen+1);
            for(int iGenM=0; iGenM< NbinM; iGenM++){
                genM12Beg = mBeg + binWidth * iGenM;
                genM12End = mBeg + binWidth * (iGenM+1);
                xM13[iGen*NbinM+iGenM] = genM13Beg + 0.5 * binM13;
                yM12[iGen*NbinM+iGenM] = genM12Beg + 0.5 * binWidth;
            }
        }


        /**
         * 1. Get Entries at generator level
         * */
        tc->SetBranchAddress("hel1213",     &costheta);
        tc->SetBranchAddress("m12",     &m12);
        cout<<" read generator level truth "<<endl;
        for(int itc=0; itc< tc->GetEntries(); itc++){
            tc->GetEntry(itc);
            if(costheta < cosBeg || costheta > cosEnd) continue;
            if(m12 < mBeg || m12 > mEnd) continue;
            m13Idx =( costheta - cosBeg)/binCos;
            m12Idx =( m12 - mBeg)/binWidth;
            NGen[m13Idx*NbinM+m12Idx]++;
        }
        cout<<" "<<__FILE__<<" "<<__LINE__<<endl;

        /**
         * 2. Get Entries after selection and efficiency matrix
         * */
        TH1F *hGen=new TH1F("hGen","hGen",100, mBeg, mEnd);
        //TFile *f =new TFile("../alg_revised/ChooseBestSigTagPair/root/0001_abs.root");
        TFile *fM =new TFile(fPhsp);
        TFile *tem=new TFile("../sample_matrix_result/tmp.root","recreate"); 
        cout<<" "<<__FILE__<<" "<<__LINE__<<endl;
        TTree* t0=(TTree*) fM->Get("DsDecay");
        cout<<" "<<__FILE__<<" "<<__LINE__<<endl;
        sprintf(cut,"isSig==1&&tagMatch==1&&m12 > %f && m12 <%f",mBeg, mEnd);
        TTree* t=t0->CopyTree(cut);
        cout<<" "<<__FILE__<<" "<<__LINE__<<endl;
        //TTree* t=(TTree*) fM->Get("KKPi");
        int *Nij= new int[NbinCos*NbinM];
        cout<<" "<<__FILE__<<" "<<__LINE__<<endl;
        int Nent= t->Project("hGen","m12",cut);
        cout<<" "<<__FILE__<<" "<<__LINE__<<endl;
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
            //if(it%1000==0) cout<<" "<<it<<" events processed"<<endl;
            m13Idx =( costheta - cosBeg)/binCos;
            m13IdxTruth = ( costheta_truth - cosBeg)/binCos;
            m12Idx =( m12 - mBeg)/binWidth;
            m12IdxTruth = ( m12_truth - mBeg)/binWidth;
            if(m13Idx<0|| m13Idx>=NbinCos || m12Idx<0 || m12Idx>=NbinM){
                cout<<" costheta = "<<costheta<<" m12 = "<<m12<<endl;
                continue;
            }
            Ni[m13Idx*NbinM+m12Idx][m13IdxTruth*NbinM + m12IdxTruth]++;
        }
        double *eff = new double[NbinCos * NbinM];
        ofstream ofs0("../sample_matrix_result/eff0array.dat");
        for(int rowIdx = 0; rowIdx < NbinCos*NbinM; rowIdx++){
            for(int rowIdxM=0; rowIdxM < NbinCos*NbinM; rowIdxM++){
                //cout<<__LINE__<<" Ni["<<rowIdx<<"]["<<rowIdxM<<"]= "<<Ni[rowIdx][rowIdxM]<<" NGen["<<rowIdxM<<"]="<<NGen[rowIdxM]<<endl;
                eff0Array[rowIdx*NbinM*NbinCos+rowIdxM]= Ni[rowIdx][rowIdxM]*1.0/NGen[rowIdxM];
                ofs0<<eff0Array[rowIdx*NbinM*NbinCos+rowIdxM]<<endl;
                if(isGenerator){
                    cout<<"isGenerator: "<<isGenerator<<" Ni["<<rowIdx<<"]= "<<Ni[rowIdx]<<" cut: "<<cut<<endl;
                    if( rowIdx == rowIdxM){
                        eff0Array[rowIdx*NbinM*NbinCos+rowIdxM]=1;
                    }
                    else{
                        eff0Array[rowIdx*NbinM*NbinCos+rowIdxM]=0;
                    }
                }
            }
        }
        ofs0.close();
        /*delete tc;
        delete t;
        delete[] Nij;
        delete[] eff;
        delete tem;
        */
    }
    else{
        ifstream ifs("../sample_matrix_result/eff0array.dat");
        for(int ni=0; ni < NbinM*NbinCos*NbinM*NbinCos; ni++){
            ifs >>  eff0Array[ni];
        }
    }
    /**
     * 3. invert the eff matrix
     * */
    cout<<" calculate inverse of eff matrix "<<endl;
    TMatrixD * mEff= new TMatrixD(NbinM*NbinCos, NbinM*NbinCos, eff0Array);
    TDecompLU *lu=new TDecompLU(*mEff);
    inverseMEff=new TMatrixD(lu->Invert());

    /**
     * 4. get efficiency array
     * */
    double *NDSel =new double[NbinM*NbinCos];
    std::uninitialized_fill_n(NDSel, NbinM*NbinCos, 0);
    cout<<" read date after selection "<<endl;
    TChain *tD0 =new TChain("DsDecay","DsDecay");
    tD0->AddFile(fDName);
    //sprintf(cut,"m12 > %f && m12 <%f",mBeg, mEnd);
    sprintf(cut,"isSig==1&&tagMatch==1&&m12 > %f && m12 <%f",mBeg, mEnd);
    TTree* tD=tD0->CopyTree(cut);
    tD->SetBranchAddress("hel1213",     &costheta);
    tD->SetBranchAddress("m12",     &m12);
    //tD->SetBranchAddress("hel1213_truth",     &costheta_truth);
    //tD->SetBranchAddress("m12_truth",     &m12_truth);
    for(int itd=0; itd< tD->GetEntries(); itd++){
        tD->GetEntry(itd);
        if(costheta < cosBeg || costheta > cosEnd) continue;
        if(m12 < mBeg || m12 > mEnd) continue;
        m13Idx =( costheta - cosBeg)/binCos;
        m12Idx =( m12 - mBeg)/binWidth;
        NDSel[m13Idx*NbinM+m12Idx]++;
    }
    cout<<" calculate eff array "<<endl;
    TMatrixD *mDSel=new TMatrixD(NbinM*NbinCos, 1, NDSel);
    TMatrixD *mDGen=new TMatrixD((*inverseMEff)*(*mDSel));
    double *inEff = inverseMEff->GetMatrixArray();
    double *aDSel=mDSel->GetMatrixArray(); 
    double *aDGen=mDGen->GetMatrixArray(); 
    double *effArray = new double[NbinM*NbinCos];
    for( int nB=0; nB < NbinM*NbinCos; nB++){
        if(fabs(aDSel[nB])< 1e-30){
            effArray[nB] = 0;
        }
        else{
            effArray[nB] = aDGen[nB]/aDSel[nB];
        }
        //cout<<" effArray["<<nB<<"]="<<effArray[nB]<<" aDSel["<<nB<<"]="<<aDSel[nB]<<" aDGen["<<nB<<"]="<<aDGen[nB]<<endl;
        /*cout<<" row: "<<nB;
        for( int nI=0; nI < NbinM*NbinCos; nI++){
            cout<<" "<<inEff[nB * NbinM*NbinCos+nI];
        }
        cout<<endl;
        */
    }
    ofstream ofs("../sample_matrix_result/effArray.dat");
    for( int nB=0; nB < NbinM*NbinCos; nB++){
        //cout<<" effArray["<<nB<<"]="<<effArray[nB]<<" aDSel["<<nB<<"]="<<aDSel[nB]<<" aDGen["<<nB<<"]="<<aDGen[nB]<<endl;
        ofs << effArray[nB]<<endl;
    }
    //MEff->Invert(effArray);
    cout<<endl;
    ofs.close();
    cout<<" start to  calculate matrix"<<endl;
    return effArray;
    //return eff0Array ;
}

void correctData(int recal=0){
    char cut[500];
    binWidth = (mEnd -mBeg)/NbinM;
    binM13 = (m13End -m13Beg)/NbinM;
    binCos = (cosEnd -cosBeg)/NbinCos;


    double *eff0Array=new double[NbinM*NbinCos*NbinM*NbinCos];
    std::uninitialized_fill_n(eff0Array, NbinM*NbinCos*NbinM*NbinCos, 0);
    double **Ni  =new double*[NbinM*NbinCos];
    for(int ii =0; ii < NbinM*NbinCos; ii++){
        Ni[ii] = new double[NbinM*NbinCos];
        std::uninitialized_fill_n(Ni[ii], NbinM*NbinCos, 0);
    }
    //int *Ni =new int[NbinM*NbinCos];
    //std::uninitialized_fill_n(Ni, NbinM*NbinCos, 0);
    double *NGen =new double[NbinM*NbinCos];
    std::uninitialized_fill_n(NGen, NbinM*NbinCos, 0);
    double *xM13 =new double[NbinM*NbinCos];
    double *yM12 =new double[NbinM*NbinCos];
    TMatrixD *inverseMEff;
    double *effArray=calEff(eff0Array, Ni, NGen, xM13, yM12, inverseMEff);
    /*ofstream ofs("../sample_matrix_../sample_matrix_result/array.dat");
      ofs<<"\t xM13\t\t yM12\t\t Ni\t\t NGen\t\t eff0Array\t"<<endl;
      for(int i =0; i< NbinM*NbinCos; i++){
      ofs<<"\t "<<xM13[i]<<"\t "<<yM12[i]<<"\t\t "<<NGen[i]<<"\t"<<effArray[i]<<endl;
    //ofs<<"\t "<<xM13[i]<<"\t "<<yM12[i]<<"\t\t "<<Ni[i]<<"\t "<<NGen[i]<<"\t"<<effArray[i]<<endl;
    }
    ofs.close();
    */
}


int main(int argc, char **argv){
    if(argc > 1){
        correctData(atoi(argv[1]));
    }
    else{
        correctData();
    }
}
