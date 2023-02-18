#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int TMVAClassification( TString myMethodList = "" )
{
   
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

      Use["BDTG"]            = 1; // uses Gradient Boost
   		Use["BDT"]            = 1;
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
 //  Use["RuleFit"]         = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   

   // ------------------------------------------------------------------------------------------------------------

   // Here the preparation phase begins

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   TChain *ch1 = new TChain("tree");
   ch1->Add("mccutsf.root");
   TTree *tr1 = ch1;

   TChain *ch2 = new TChain("tree");
   ch2->Add("datacutsf_r.root");
   TTree *tr2 = ch2;

   // Register the training and test trees

   TTree *signalTree     = tr1;
   TTree *background     = tr2;

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE");

   
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
  
   dataloader->AddVariable( "maxkpt:=TMath::Max(Kmpt, Kppt)",  'F' );
   dataloader->AddVariable( "maxktrk:=TMath::Max(KmtrkMinIP/KmtrkMinIPE, KptrkMinIP/KptrkMinIPE)", 'F' );
  dataloader->AddVariable( "maxmu:=TMath::Max(MumMinIP/MumMinIPE, MupMinIP/MupMinIPE)", 'F' );
   dataloader->AddVariable( "Blxysig", 'F' );
   dataloader->AddVariable( "Bvtxcl", 'F' );
//  dataloader->AddVariable( "maxkmtrkdca:=TMath::Max(Kmtrkdcasigbs, Kptrkdcasigbs)", 'F' );
   dataloader->AddVariable( "Bcosalphabs2d", 'F' );
    dataloader->AddVariable( "Bsdcasigbs", 'F' );
	dataloader->AddVariable( "maxktriso:=TMath::Max(kmtrkIso, kptrkIso)", 'F' );  


   
   
   

   

   dataloader->AddSpectator( "Bmass",  "Spectator 1", "units", 'F' );
   dataloader->AddSpectator( "Mumumass", "Spectator 2", "units", 'F' );


   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 35.9/2788.4;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( background, backgroundWeight );

   
   // Apply additional cuts on the signal and background samples (can be different)
   //dataloader->SetBackgroundWeightExpression( "weight" );
  TCut mycuts="Blxysig>0 && Blxysig<9999 && Bsdcasigbs>0 && Bsdcasigbs<9999 && TMath::Max(MumMinIP/MumMinIPE, MupMinIP/MupMinIPE)>0 && TMath::Max(KmtrkMinIP/KmtrkMinIPE, KptrkMinIP/KptrkMinIPE)>0 ";//&& TMath::Max(MumMinIP/MumMinIPE, MupMinIP/MupMinIPE)<9999 && TMath::Max(KmtrkMinIP/KmtrkMinIPE, KptrkMinIP/KptrkMinIPE)>9999";
   
   TCut mycutb="Blxysig>0 && Blxysig<9999 && Bsdcasigbs>0 && Bsdcasigbs<9999   && TMath::Max(MumMinIP/MumMinIPE, MupMinIP/MupMinIPE)>0 && TMath::Max(KmtrkMinIP/KmtrkMinIPE, KptrkMinIP/KptrkMinIPE)>0 ";//&& TMath::Max(MumMinIP/MumMinIPE, MupMinIP/MupMinIPE)<9999 && TMath::Max(KmtrkMinIP/KmtrkMinIPE, KptrkMinIP/KptrkMinIPE)>9999"; // for example: TCut mycutb = "abs(var1)<0.5";

   
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   

   // Cut optimisation
   
   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   
  
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   return TMVAClassification(methodList);
}
