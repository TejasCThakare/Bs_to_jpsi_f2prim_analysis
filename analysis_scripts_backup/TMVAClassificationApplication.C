/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained classifiers
/// within an analysis module
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Exectuable: TMVAClassificationApplication
///
//   run:
//   root -l TMVAClassificationApplication.C\(\"BDT\"\)
//
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void TMVAClassificationApplication( TString myMethodList = "" )
{

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   Use["DNN_CPU"]         = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"]         = 0;         // Multi-core accelerated DNN.
   //
   // Support Vector Machine
   Use["SVM"]             = 0;
   //
   // Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t var1, var2, var3, var4, var5, var6, var7, var8 ;
   Int_t IntVar1;

   


   reader->AddVariable( "maxkpt:=TMath::Max(Kmpt, Kppt)", &var1);
   reader->AddVariable( "maxktrk:=TMath::Max(KmtrkMinIP/KmtrkMinIPE, KptrkMinIP/KptrkMinIPE)",             &var2);
   reader->AddVariable( "maxmu:=TMath::Max(MumMinIP/MumMinIPE, MupMinIP/MupMinIPE)",       &var3);
   reader->AddVariable( "Blxysig", &var4);
   reader->AddVariable( "Bvtxcl",   &var5);
   reader->AddVariable( "Bcosalphabs2d",  &var6 );
   reader->AddVariable( "Bsdcasigbs", &var7);
   reader->AddVariable( "maxktriso:=TMath::Max(kmtrkIso, kptrkIso)",   &var8 );
  // reader->AddVariable( "Bsdcasigbs",   &var9);
  // reader->AddVariable( "Phimass", &var10);
  // reader->AddVariable( "BsIso",   &var11);
  // reader->AddVariable( "K_Iso := TMath::Max(kmtrkIso, kptrkIso)",  &var12);

   Float_t spec1,spec2;//spec3;                                                                                            
   reader->AddSpectator( "Bmass",           &spec1 );                                         
   reader->AddSpectator( "Mumumass",        &spec2 ); 
   //reader->AddSpectator( "Mumumasserr",     &spec3 ); 


   //  Book the MVA methods

   TString dir    = "dataset/weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 100;
   TH1F *histLk(0);
   TH1F *histLkD(0);
   TH1F *histLkPCA(0);
   TH1F *histLkKDE(0);
   TH1F *histLkMIX(0);
   TH1F *histPD(0);
   TH1F *histPDD(0);
   TH1F *histPDPCA(0);
   TH1F *histPDEFoam(0);
   TH1F *histPDEFoamErr(0);
   TH1F *histPDEFoamSig(0);
   TH1F *histKNN(0);
   TH1F *histHm(0);
   TH1F *histFi(0);
   TH1F *histFiG(0);
   TH1F *histFiB(0);
   TH1F *histLD(0);
   TH1F *histNn(0);
   TH1F *histNnbfgs(0);
   TH1F *histNnbnn(0);
   TH1F *histNnC(0);
   TH1F *histNnT(0);
   TH1F *histBdt(0);
   TH1F *histBdtG(0);
   TH1F *histBdtB(0);
   TH1F *histBdtD(0);
   TH1F *histBdtF(0);
   TH1F *histRf(0);
   TH1F *histSVMG(0);
   TH1F *histSVMP(0);
   TH1F *histSVML(0);
   TH1F *histFDAMT(0);
   TH1F *histFDAGA(0);
   TH1F *histCat(0);
   TH1F *histPBdt(0);
   TH1F *histDnnGpu(0);
   TH1F *histDnnCpu(0);

  
  

   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );

  

   
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.



   //****************************                                                                                                                    
   //  define input trees                                                                                                                             
   //****************************                                                                                                                     
   TChain *ch_sig = new TChain("tree");
   ch_sig->Add("sel_BsToJpsif2p_2016MC_UL_Official_Presel_mc.lite_cutopt.root");
   TTree *sigtr_in = ch_sig;
   std::cout << " @@@  entries in sig tree: " << sigtr_in->GetEntries() << std::endl;

   TChain *ch_bkg = new TChain("tree");
   ch_bkg->Add("sel_Combine_2016_Mini_Presel_data_cutopt.root");
   TTree *bkgtr_in = ch_bkg;
   std::cout << " @@@  entries in bkg tree: " << bkgtr_in->GetEntries() << std::endl;



   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //


   std::cout << "--- Select signal sample" << std::endl;

   Double_t userIntVar4, userVar1, userVar2, userVar3, userVar4, userVar5, userVar6, userVar7, userVar8, userVar9, userVar10, userVar11, userVar12, userVar13, userVar14, userVar15, userVar16, userVar17, userVar18, userVar19, userVar20, userVar21, userVar22, userVar23, userVar24, userVar25, userVar26, userVar27, userVar28, userVar29, userVar30,  userVar31, userVar32, userVar33, userVar34, userVar35, userVar36, userVar37, userVar38, userVar39, userVar40, userVar41, userVar42, userVar43, userVar44, userVar45, userVar46, userVar47, userVar48, userVar49, userVar50, userVar51, userVar52, userVar53, userVar54, userVar55, userVar56, userVar57, userVar60 , userVar61 ;
   Int_t userIntVar1, userIntVar2, userIntVar3, userVar58, userVar59, userVar62;

   sigtr_in->SetBranchAddress( "Mumumass",          &userVar1 );
   sigtr_in->SetBranchAddress( "Mumumasserr",          &userVar2 );
   sigtr_in->SetBranchAddress( "Phimass", &userVar3 );
   sigtr_in->SetBranchAddress( "Kmpt", &userVar4 );
   sigtr_in->SetBranchAddress( "Kppt",        &userVar5 );
   sigtr_in->SetBranchAddress( "Kmeta",       &userVar6 );
   sigtr_in->SetBranchAddress( "Kpeta", &userVar7 );
   sigtr_in->SetBranchAddress( "Kmphi", &userVar8 );
   sigtr_in->SetBranchAddress( "Kpphi",         &userVar9 );
   sigtr_in->SetBranchAddress( "Kmtrkdcasigbs",      &userVar10);
   sigtr_in->SetBranchAddress( "Kptrkdcasigbs",   &userVar11);
   sigtr_in->SetBranchAddress( "Mumpt",       &userVar12);   
   sigtr_in->SetBranchAddress( "Muppt",           &userVar13);  
   sigtr_in->SetBranchAddress( "Mumeta",           &userVar14);
   sigtr_in->SetBranchAddress( "Mupeta",    &userVar15);
   sigtr_in->SetBranchAddress( "Mumphi",    &userVar16);
   sigtr_in->SetBranchAddress( "Mupphi",      &userVar17);
   sigtr_in->SetBranchAddress( "Mumdcasigbs",      &userVar18);
   sigtr_in->SetBranchAddress( "Mupdcasigbs",    &userVar19);
   sigtr_in->SetBranchAddress( "Bsdcasigbs",         &userVar20);
   sigtr_in->SetBranchAddress( "MumMinIP",      &userVar21);
   sigtr_in->SetBranchAddress( "MupMinIP",      &userVar22);
   sigtr_in->SetBranchAddress( "MumMinIPE",     &userVar23);
   sigtr_in->SetBranchAddress( "MupMinIPE",     &userVar24);
   sigtr_in->SetBranchAddress( "KptrkMinIP",   &userVar25);
   sigtr_in->SetBranchAddress( "KptrkMinIPE",   &userVar26);
   sigtr_in->SetBranchAddress( "KmtrkMinIP",   &userVar27);
   sigtr_in->SetBranchAddress( "KmtrkMinIPE",   &userVar28);
   sigtr_in->SetBranchAddress( "MumMinIP2D",           &userVar29);
   sigtr_in->SetBranchAddress( "MupMinIP2D",    &userVar30);
   sigtr_in->SetBranchAddress( "MumMinIP2DE",    &userVar31);
   sigtr_in->SetBranchAddress( "MupMinIP2DE",    &userVar32);
   sigtr_in->SetBranchAddress( "mupIso", &userVar33 );
   sigtr_in->SetBranchAddress( "mumIso", &userVar34 );
   sigtr_in->SetBranchAddress( "BsIso",        &userVar35 );
   sigtr_in->SetBranchAddress( "kptrkIso",       &userVar36 );
   sigtr_in->SetBranchAddress( "kmtrkIso", &userVar37 );
   sigtr_in->SetBranchAddress( "Bmass", &userVar38 );
   sigtr_in->SetBranchAddress( "Bpt",         &userVar39 );
   sigtr_in->SetBranchAddress( "Beta",      &userVar40);
   sigtr_in->SetBranchAddress( "Bphi",   &userVar41);
   sigtr_in->SetBranchAddress( "Phipt",       &userVar42);   
   sigtr_in->SetBranchAddress( "Phieta",           &userVar43);  
   sigtr_in->SetBranchAddress( "Phiphi",           &userVar44);
   sigtr_in->SetBranchAddress( "Bvtxcl",    &userVar45);
   sigtr_in->SetBranchAddress( "Blxysig",    &userVar46);
   sigtr_in->SetBranchAddress( "Bcosalphabs",      &userVar47);
   sigtr_in->SetBranchAddress( "Bcosalphabs2d",      &userVar48);
   sigtr_in->SetBranchAddress( "Q2",    &userVar49);
   sigtr_in->SetBranchAddress( "dimupt",         &userVar50);
   sigtr_in->SetBranchAddress( "dimueta",      &userVar51);
   sigtr_in->SetBranchAddress( "dimuvtxcl",      &userVar52);
   sigtr_in->SetBranchAddress( "dimulsig",     &userVar53);
   sigtr_in->SetBranchAddress( "dimuDCA",     &userVar54);
   sigtr_in->SetBranchAddress( "CosThetaL",   &userVar55);
   sigtr_in->SetBranchAddress( "CosThetaK",   &userVar56);
   sigtr_in->SetBranchAddress( "Phi",           &userVar57);
   sigtr_in->SetBranchAddress( "ptrkqual",    &userVar58);
   sigtr_in->SetBranchAddress( "mtrkqual",    &userVar59);
   sigtr_in->SetBranchAddress( "dr0",    &userVar60);
   sigtr_in->SetBranchAddress( "dr1",           &userVar61);
 //  sigtr_in->SetBranchAddress( "eventID",           &userVar62);   
//
   sigtr_in->SetBranchAddress( "JpsiTriggers",  &userIntVar1); 
   sigtr_in->SetBranchAddress( "PsiPTriggers",  &userIntVar2);
   sigtr_in->SetBranchAddress( "LMNTTriggers",  &userIntVar3);
  // sigtr_in->SetBranchAddress( "Puw8",      &userIntVar4);
   
    
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   bkgtr_in->SetBranchAddress( "Mumumass",          &userVar1 );
   bkgtr_in->SetBranchAddress( "Mumumasserr",          &userVar2 );
   bkgtr_in->SetBranchAddress( "Phimass", &userVar3 );
   bkgtr_in->SetBranchAddress( "Kmpt", &userVar4 );
   bkgtr_in->SetBranchAddress( "Kppt",        &userVar5 );
   bkgtr_in->SetBranchAddress( "Kmeta",       &userVar6 );
   bkgtr_in->SetBranchAddress( "Kpeta", &userVar7 );
   bkgtr_in->SetBranchAddress( "Kmphi", &userVar8 );
   bkgtr_in->SetBranchAddress( "Kpphi",         &userVar9 );
   bkgtr_in->SetBranchAddress( "Kmtrkdcasigbs",      &userVar10);
   bkgtr_in->SetBranchAddress( "Kptrkdcasigbs",   &userVar11);
   bkgtr_in->SetBranchAddress( "Mumpt",       &userVar12);   
   bkgtr_in->SetBranchAddress( "Muppt",           &userVar13);  
   bkgtr_in->SetBranchAddress( "Mumeta",           &userVar14);
   bkgtr_in->SetBranchAddress( "Mupeta",    &userVar15);
   bkgtr_in->SetBranchAddress( "Mumphi",    &userVar16);
   bkgtr_in->SetBranchAddress( "Mupphi",      &userVar17);
   bkgtr_in->SetBranchAddress( "Mumdcasigbs",      &userVar18);
   bkgtr_in->SetBranchAddress( "Mupdcasigbs",    &userVar19);
   bkgtr_in->SetBranchAddress( "Bsdcasigbs",         &userVar20);
   bkgtr_in->SetBranchAddress( "MumMinIP",      &userVar21);
   bkgtr_in->SetBranchAddress( "MupMinIP",      &userVar22);
   bkgtr_in->SetBranchAddress( "MumMinIPE",     &userVar23);
   bkgtr_in->SetBranchAddress( "MupMinIPE",     &userVar24);
   bkgtr_in->SetBranchAddress( "KptrkMinIP",   &userVar25);
   bkgtr_in->SetBranchAddress( "KptrkMinIPE",   &userVar26);
   bkgtr_in->SetBranchAddress( "KmtrkMinIP",   &userVar27);
   bkgtr_in->SetBranchAddress( "KmtrkMinIPE",   &userVar28);
   bkgtr_in->SetBranchAddress( "MumMinIP2D",           &userVar29);
   bkgtr_in->SetBranchAddress( "MupMinIP2D",    &userVar30);
   bkgtr_in->SetBranchAddress( "MumMinIP2DE",    &userVar31);
   bkgtr_in->SetBranchAddress( "MupMinIP2DE",    &userVar32);
   bkgtr_in->SetBranchAddress( "mupIso", &userVar33 );
   bkgtr_in->SetBranchAddress( "mumIso", &userVar34 );
   bkgtr_in->SetBranchAddress( "BsIso",        &userVar35 );
   bkgtr_in->SetBranchAddress( "kptrkIso",       &userVar36 );
   bkgtr_in->SetBranchAddress( "kmtrkIso", &userVar37 );
   bkgtr_in->SetBranchAddress( "Bmass", &userVar38 );
   bkgtr_in->SetBranchAddress( "Bpt",         &userVar39 );
   bkgtr_in->SetBranchAddress( "Beta",      &userVar40);
   bkgtr_in->SetBranchAddress( "Bphi",   &userVar41);
   bkgtr_in->SetBranchAddress( "Phipt",       &userVar42);   
   bkgtr_in->SetBranchAddress( "Phieta",           &userVar43);  
   bkgtr_in->SetBranchAddress( "Phiphi",           &userVar44);
   bkgtr_in->SetBranchAddress( "Bvtxcl",    &userVar45);
   bkgtr_in->SetBranchAddress( "Blxysig",    &userVar46);
   bkgtr_in->SetBranchAddress( "Bcosalphabs",      &userVar47);
   bkgtr_in->SetBranchAddress( "Bcosalphabs2d",      &userVar48);
   bkgtr_in->SetBranchAddress( "Q2",    &userVar49);
   bkgtr_in->SetBranchAddress( "dimupt",         &userVar50);
   bkgtr_in->SetBranchAddress( "dimueta",      &userVar51);
   bkgtr_in->SetBranchAddress( "dimuvtxcl",      &userVar52);
   bkgtr_in->SetBranchAddress( "dimulsig",     &userVar53);
   bkgtr_in->SetBranchAddress( "dimuDCA",     &userVar54);
   bkgtr_in->SetBranchAddress( "CosThetaL",   &userVar55);
   bkgtr_in->SetBranchAddress( "CosThetaK",   &userVar56);
   bkgtr_in->SetBranchAddress( "Phi",           &userVar57);
   bkgtr_in->SetBranchAddress( "ptrkqual",    &userVar58);
   bkgtr_in->SetBranchAddress( "mtrkqual",    &userVar59);
   bkgtr_in->SetBranchAddress( "dr0",    &userVar60);
   bkgtr_in->SetBranchAddress( "dr1",           &userVar61);
  // bkgtr_in->SetBranchAddress( "eventID",           &userVar62);
//
   bkgtr_in->SetBranchAddress( "JpsiTriggers",  &userIntVar1); 
   bkgtr_in->SetBranchAddress( "PsiPTriggers",  &userIntVar2);
   bkgtr_in->SetBranchAddress( "LMNTTriggers",  &userIntVar3);
  // bkgtr_in->SetBranchAddress( "Puw8",      &userIntVar4);
  

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;


   TFile *targetS  = new TFile( "sel_BsToJpsif2p_2016MC_UL_Official_Presel_mc.lite_cutopt_bdtr.root","RECREATE" );
   TTree *treeS    = new TTree("tree","sig tree");
//------------treeS->Branch()
   treeS->Branch( "Mumumass",          &userVar1 );
   treeS->Branch( "Mumumasserr",          &userVar2 );
   treeS->Branch( "Phimass", &userVar3 );
   treeS->Branch( "Kmpt", &userVar4 );
   treeS->Branch( "Kppt",        &userVar5 );
   treeS->Branch( "Kmeta",       &userVar6 );
   treeS->Branch( "Kpeta", &userVar7 );
   treeS->Branch( "Kmphi", &userVar8 );
   treeS->Branch( "Kpphi",         &userVar9 );
   treeS->Branch( "Kmtrkdcasigbs",      &userVar10);
   treeS->Branch( "Kptrkdcasigbs",   &userVar11);
   treeS->Branch( "Mumpt",       &userVar12);   
   treeS->Branch( "Muppt",           &userVar13);  
   treeS->Branch( "Mumeta",           &userVar14);
   treeS->Branch( "Mupeta",    &userVar15);
   treeS->Branch( "Mumphi",    &userVar16);
   treeS->Branch( "Mupphi",      &userVar17);
   treeS->Branch( "Mumdcasigbs",      &userVar18);
   treeS->Branch( "Mupdcasigbs",    &userVar19);
   treeS->Branch( "Bsdcasigbs",         &userVar20);
   treeS->Branch( "MumMinIP",      &userVar21);
   treeS->Branch( "MupMinIP",      &userVar22);
   treeS->Branch( "MumMinIPE",     &userVar23);
   treeS->Branch( "MupMinIPE",     &userVar24);
   treeS->Branch( "KptrkMinIP",   &userVar25);
   treeS->Branch( "KptrkMinIPE",   &userVar26);
   treeS->Branch( "KmtrkMinIP",   &userVar27);
   treeS->Branch( "KmtrkMinIPE",   &userVar28);
   treeS->Branch( "MumMinIP2D",           &userVar29);
   treeS->Branch( "MupMinIP2D",    &userVar30);
   treeS->Branch( "MumMinIP2DE",    &userVar31);
   treeS->Branch( "MupMinIP2DE",    &userVar32);
   treeS->Branch( "mupIso", &userVar33 );
   treeS->Branch( "mumIso", &userVar34 );
   treeS->Branch( "BsIso",        &userVar35 );
   treeS->Branch( "kptrkIso",       &userVar36 );
   treeS->Branch( "kmtrkIso", &userVar37 );
   treeS->Branch( "Bmass", &userVar38 );
   treeS->Branch( "Bpt",         &userVar39 );
   treeS->Branch( "Beta",      &userVar40);
   treeS->Branch( "Bphi",   &userVar41);
   treeS->Branch( "Phipt",       &userVar42);   
   treeS->Branch( "Phieta",           &userVar43);  
   treeS->Branch( "Phiphi",           &userVar44);
   treeS->Branch( "Bvtxcl",    &userVar45);
   treeS->Branch( "Blxysig",    &userVar46);
   treeS->Branch( "Bcosalphabs",      &userVar47);
   treeS->Branch( "Bcosalphabs2d",      &userVar48);
   treeS->Branch( "Q2",    &userVar49);
   treeS->Branch( "dimupt",         &userVar50);
   treeS->Branch( "dimueta",      &userVar51);
   treeS->Branch( "dimuvtxcl",      &userVar52);
   treeS->Branch( "dimulsig",     &userVar53);
   treeS->Branch( "dimuDCA",     &userVar54);
   treeS->Branch( "CosThetaL",   &userVar55);
   treeS->Branch( "CosThetaK",   &userVar56);
   treeS->Branch( "Phi",           &userVar57);
   treeS->Branch( "ptrkqual",    &userVar58);
   treeS->Branch( "mtrkqual",    &userVar59);
   treeS->Branch( "dr0",    &userVar60);
   treeS->Branch( "dr1",           &userVar61);
  // treeS->Branch( "eventID",           &userVar62);
//
   treeS->Branch( "JpsiTriggers",  &userIntVar1); 
   treeS->Branch( "PsiPTriggers",  &userIntVar2);
   treeS->Branch( "LMNTTriggers",  &userIntVar3);
  // treeS->Branch( "Puw8",      &userIntVar4);
     
//---------------------treeS->Branch ends--------------
  
   double BDT_sig, BDTG_sig;
   treeS->Branch("BDT", &BDT_sig);
   treeS->Branch("BDTG", &BDTG_sig);


   std::cout << "--- Processing: " << sigtr_in->GetEntries() << " signal events" << std::endl;

   for (Long64_t ievt=0; ievt<sigtr_in->GetEntries(); ievt++) {

      //if (ievt==0) continue;
      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      if (ievt==1000) continue;
      
      sigtr_in->GetEntry(ievt);
      
      var1  = TMath::Max(userVar4, userVar5);
      var2  = TMath::Max(userVar27/userVar28, userVar25/userVar26);
      var3  = TMath::Max(userVar21/userVar23, userVar22/userVar24);
      var4  = userVar46;
      var5  = userVar45;
      var6  = userVar48;
      var7  = userVar20;
      var8  = TMath::Max(userVar37, userVar36);
      spec1 = userVar38;
      spec2 = userVar1;
   
     
      if (TMath::IsNaN(var1)) continue;
      if (TMath::IsNaN(var2)) continue;
      if (TMath::IsNaN(var3)) continue;
      if (TMath::IsNaN(var4)) continue;
      if (TMath::IsNaN(var5)) continue;
      if (TMath::IsNaN(var6)) continue;
      if (TMath::IsNaN(var7)) continue;
      if (TMath::IsNaN(var8)) continue;
      if (TMath::IsNaN(spec1)) continue;
	  if (TMath::IsNaN(spec2)) continue;      
 	      
      BDT_sig  = reader->EvaluateMVA( "BDT method");
      BDTG_sig = reader->EvaluateMVA( "BDTG method");

      treeS->Fill();

      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
   }

  
   std::cout << "--- End of signal event loop: " << std::endl; 

   treeS->Write();


   if (Use["BDT"          ])   histBdt    ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write();



   targetS->Close();

   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;




   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   TFile *targetB  = new TFile( "sel_Combine_2016_Mini_Presel_data_cutopt_bdtr.root","RECREATE" );
   TTree *treeB    = new TTree("tree","bkg tree");

//---------------------treeB->Branch starts-----------------
   treeB->Branch( "Mumumass",          &userVar1 );
   treeB->Branch( "Mumumasserr",          &userVar2 );
   treeB->Branch( "Phimass", &userVar3 );
   treeB->Branch( "Kmpt", &userVar4 );
   treeB->Branch( "Kppt",        &userVar5 );
   treeB->Branch( "Kmeta",       &userVar6 );
   treeB->Branch( "Kpeta", &userVar7 );
   treeB->Branch( "Kmphi", &userVar8 );
   treeB->Branch( "Kpphi",         &userVar9 );
   treeB->Branch( "Kmtrkdcasigbs",      &userVar10);
   treeB->Branch( "Kptrkdcasigbs",   &userVar11);
   treeB->Branch( "Mumpt",       &userVar12);   
   treeB->Branch( "Muppt",           &userVar13);  
   treeB->Branch( "Mumeta",           &userVar14);
   treeB->Branch( "Mupeta",    &userVar15);
   treeB->Branch( "Mumphi",    &userVar16);
   treeB->Branch( "Mupphi",      &userVar17);
   treeB->Branch( "Mumdcasigbs",      &userVar18);
   treeB->Branch( "Mupdcasigbs",    &userVar19);
   treeB->Branch( "Bsdcasigbs",         &userVar20);
   treeB->Branch( "MumMinIP",      &userVar21);
   treeB->Branch( "MupMinIP",      &userVar22);
   treeB->Branch( "MumMinIPE",     &userVar23);
   treeB->Branch( "MupMinIPE",     &userVar24);
   treeB->Branch( "KptrkMinIP",   &userVar25);
   treeB->Branch( "KptrkMinIPE",   &userVar26);
   treeB->Branch( "KmtrkMinIP",   &userVar27);
   treeB->Branch( "KmtrkMinIPE",   &userVar28);
   treeB->Branch( "MumMinIP2D",           &userVar29);
   treeB->Branch( "MupMinIP2D",    &userVar30);
   treeB->Branch( "MumMinIP2DE",    &userVar31);
   treeB->Branch( "MupMinIP2DE",    &userVar32);
   treeB->Branch( "mupIso", &userVar33 );
   treeB->Branch( "mumIso", &userVar34 );
   treeB->Branch( "BsIso",        &userVar35 );
   treeB->Branch( "kptrkIso",       &userVar36 );
   treeB->Branch( "kmtrkIso", &userVar37 );
   treeB->Branch( "Bmass", &userVar38 );
   treeB->Branch( "Bpt",         &userVar39 );
   treeB->Branch( "Beta",      &userVar40);
   treeB->Branch( "Bphi",   &userVar41);
   treeB->Branch( "Phipt",       &userVar42);   
   treeB->Branch( "Phieta",           &userVar43);  
   treeB->Branch( "Phiphi",           &userVar44);
   treeB->Branch( "Bvtxcl",    &userVar45);
   treeB->Branch( "Blxysig",    &userVar46);
   treeB->Branch( "Bcosalphabs",      &userVar47);
   treeB->Branch( "Bcosalphabs2d",      &userVar48);
   treeB->Branch( "Q2",    &userVar49);
   treeB->Branch( "dimupt",         &userVar50);
   treeB->Branch( "dimueta",      &userVar51);
   treeB->Branch( "dimuvtxcl",      &userVar52);
   treeB->Branch( "dimulsig",     &userVar53);
   treeB->Branch( "dimuDCA",     &userVar54);
   treeB->Branch( "CosThetaL",   &userVar55);
   treeB->Branch( "CosThetaK",   &userVar56);
   treeB->Branch( "Phi",           &userVar57);
   treeB->Branch( "ptrkqual",    &userVar58);
   treeB->Branch( "mtrkqual",    &userVar59);
   treeB->Branch( "dr0",    &userVar60);
   treeB->Branch( "dr1",           &userVar61);
//   treeB->Branch( "eventID",           &userVar62);
//
   treeB->Branch( "JpsiTriggers",  &userIntVar1); 
   treeB->Branch( "PsiPTriggers",  &userIntVar2);
   treeB->Branch( "LMNTTriggers",  &userIntVar3);
  // treeB->Branch( "Puw8",      &userIntVar4);//-----------------------treeB branch ends--------------------------------
                                                          
   Float_t BDT_bkg, BDTG_bkg;
   treeB->Branch("BDT", &BDT_bkg);
   treeB->Branch("BDTG", &BDTG_bkg);


   std::cout << "--- Processing: " << bkgtr_in->GetEntries() << " background events" << std::endl;


   //////TStopwatch sw;
   //////sw.Start();


   for (Long64_t ievt=0; ievt<bkgtr_in->GetEntries(); ievt++) {
   	
   	 //if (ievt==0) continue;

     if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     if (ievt==1000) continue;

     bkgtr_in->GetEntry(ievt);

      var1  = TMath::Max(userVar4, userVar5);
      var2  = TMath::Max(userVar27/userVar28, userVar25/userVar26);
      var3  = TMath::Max(userVar21/userVar23, userVar22/userVar24);
      var4  = userVar46;
      var5  = userVar45;
      var6  = userVar48;
      var7  = userVar20;
      var8  = TMath::Max(userVar37, userVar36);
      spec1 = userVar38;
      spec2 = userVar1;
      
      
      if (TMath::IsNaN(var1)) continue;
      if (TMath::IsNaN(var2)) continue;
      if (TMath::IsNaN(var3)) continue;
      if (TMath::IsNaN(var4)) continue;
      if (TMath::IsNaN(var5)) continue;
      if (TMath::IsNaN(var6)) continue;
      if (TMath::IsNaN(var7)) continue;
      if (TMath::IsNaN(var8)) continue;
      if (TMath::IsNaN(spec1)) continue;
	  if (TMath::IsNaN(spec2)) continue;    

     BDT_bkg  = reader->EvaluateMVA( "BDT method");
     BDTG_bkg = reader->EvaluateMVA( "BDTG method");

     treeB->Fill();

     if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
     if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );


   }

   std::cout << "--- End of background event loop: " << std::endl; 

   treeB->Write();


   if (Use["BDT"          ])   histBdt    ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write();


   targetB->Close();

   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}


int main( int argc, char** argv )
{
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAClassificationApplication(methodList);
   return 0;
}


