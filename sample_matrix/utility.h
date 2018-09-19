#ifndef UTITILY_H
#define UTITILY_H
#include <memory>
#include <iostream>
#include <iomanip>
#include <TTree.h>
#include <TFile.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <TChain.h>
#include "TLorentzVector.h"
#include <TH1F.h>
#include <TMatrixD.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <Math/IFunction.h>

#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TDecompLU.h"

#include <Riostream.h>
#include "TLegend.h"
#include "TLegendEntry.h"

#include "Math/IFunction.h"
#include <cmath>
#include "TSystem.h"
#include <TStyle.h>
#include <cstring>
const static int Nbin=40;
const static double mBeg=0.988;
const static double mEnd=1.115;
const static double m13Beg=0.62;
const static double m13End=1.48;
const static double hBeg=0.988;
const static double hEnd=1.115;
const static double bkgRatio = 0.05;
const static double mK = 0.493677;
const static double cosBeg = -1;
const static double cosEnd = 1;
//const static double hBeg=0.95;
//const static double hEnd=1.15;
//const static double mBeg=8.94427190999915855e-01 + 0.08;
//const static double mEnd=1.87082869338697066 - 0.025;
const static bool debug = false;
const static bool useMatrix = false;
static bool isGenerator = false;
double binWidth=0;
double binM13=0;
double binCos=0;
using std::cout;
using std::endl;
TFile *fR;
//const char* fDName="Oth_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";
//const char* fDName="f0_0.922_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";
//const char* fDName="f0_1.02_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";
//const char* fDName="f0_0.99_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";
//const char* fDName="../../CalEff_Data/f0_0.99_test/ChooseBestSignalTagPair_SingleTag/test.root";

//const char* fDName="f0_test/ChooseBestSignalTagPair_SingleTag/test.root";

//const char* fDName="f0_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";


//const char* fDName="Phi_f0_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";

//const char* fDName="../CalEff_Data/f0_0.99_test/CompareGenSel/match.root";
//const char* fDName="Phi_test/CompareGenSel/match.root";

//const char* fDName="Phi_test/CompareGenSel/matchSel.root";
//const char* fDName="../../CalEff_Data/Phi_test/ChooseBestSignalTagPair_SingleTag/test.root";


//const char* fDName="Phi_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";

//const char* fDName="../../generic_703/GeneratorLevel_SingleTag_Recheck/rootFile/singleTag.root";
//const char* fDName="../../generic_703/GeneratorLevel_SingleTag_Recheck/rootFile/singleTag.root";
//const char* fDName="f0980_test/ana/singleTag.root";
//const char* fDName="../CalEff_Data/f0980_test/ChooseBestSignalTagPair_SingleTag/test.root";

//const char* fDName="f0980_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";
//const char* fDName="../../generic_703/alg/ChooseBestSigTagPair_SingleTag/test.root";
//const char* fDName="../root/data/test.root";
const char* fDName="../root/generic_703/test.root";
//const char* fDName="../../CalEff_Data/../../generic_703/alg/ChooseBestSigTagPair_SingleTag/test.root";
//const char* fDName="/besfs/groups/psipp/psippgroup/public/wangmeng/DsTag/KKPi/DongLY/DongLY/data_703/CalEff_Data/f0980_test/ChooseBestSignalTagPair_SingleTag/test.root";
//const char* fDName="../../data_703/alg_revised/ana/singleTag.root";
//const char* fDName="/besfs/groups/psipp/psippgroup/public/wangmeng/DsTag/KKPi/DongLY/DongLY/IO/DIY_MC_5c_cut_P7/alg_703_New_afRatio_Lagrange_BES2/ana/result.root";

//const char* fDName="/besfs/groups/psipp/psippgroup/public/wangmeng/DsTag/KKPi/DongLY/DongLY/generic_703/alg/ana/result.root";
//const char* fGen="../../generic_703/GeneratorLevel_SingleTag/rootFile/singleTag.root";
//const char* fGen="Phi_test/GeneratorLevel_SingleTag/rootFile/singleTag.root";

//const char* fGen= "Phi_test/PHSP_VER/GeneratorLevel_SingleTag/rootFile/singleTag.root";

//const char* fPhsp="Phi_test/PHSP_VER/ChooseBestSignalTagPair_SingleTag/test.root";

const char* fGen="../root/phsp/singleTag.root";
const char* fPhsp="../root/phsp/phsp_test.root";
//const char* fPhsp="../../phsp_703/alg_revised/ana/singleTag.root";
//const char* fPhsp="../../generic_703/alg/ChooseBestSigTagPair_SingleTag/test.root";
//const char* fPhsp="Phi_test/ChooseBestSignalTagPair_SingleTag/test.root";

//const char* fPhsp="../../phsp_703/alg_revised/ChooseBestSigTagPair_SingleTag/test.root";
const char* fSig="../root/generic_703/ana_singleTag.root";
//const char* fSig="../alg_revised/ana/result.root";

#endif
