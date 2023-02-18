#include <TChain.h>
#include <TChainElement.h>

void printListOfTChainElements(TChain *chain){
  TObjArray *fileElements=chain->GetListOfFiles();
  int nFiles = fileElements->GetEntries();                                                        
  TIter next(fileElements);
  TChainElement *chEl=0;
  for( int entry=0; entry < nFiles; entry++ ) {
    chEl=(TChainElement*)next();
    printf("%s\n",chEl->GetTitle());
  }
  printf("DEBUG\t\t: %d files in the chain\n",nFiles);
}


void selectSIGNAL_LMNR(const char* outfile, const char* path, const char* file, TCut cuts /*,const int flag*/){ 

                                                                             
  const std::string treename  = "tree";
  TCut mycuts  =  cuts;
  
  TChain *ch_sig = new TChain("tree");
  ch_sig->Add(Form("%s/%s", path, file));
  printListOfTChainElements(ch_sig);
  TTree *inTree = ch_sig;

  if( !inTree ) std::cout << "tree " << treename << " does not exist" << std::endl;
  std::cout << " @@@  entries in sig tree: " << inTree->GetEntries() << std::endl;

 
      
double n_pre = inTree->GetEntries();
double n_pre_wcut = inTree->GetEntries(mycuts);

std::cout << "# of events in original tuple = " << n_pre << endl;  
std::cout << "# of events in original tuple after cuts = " << n_pre_wcut << endl;  

std::cout << "Copying tree: " << treename.c_str() << endl;
    
TFile* newFile = new TFile(outfile,"RECREATE");
TTree* outTree = inTree->CopyTree( mycuts );
    
int percentCounter = 1;

for(int i = 0; i < inTree->GetEntries(); ++i){
  
  const int percent = (int)(inTree->GetEntries()/100.0);
    
  if( i == percent*percentCounter ){
    std::cout << percentCounter << " %" << std::endl;
    percentCounter++;
  }
      
 }
    
outTree->Write();
newFile->Save();

double n_post = outTree->GetEntries();
    
std::cout << "# of events in the signal tuple = " << n_post << endl;


} 

void readbranch()
{
    const char* path="/eos/user/t/tthakare/projectWork/exe_files";
    
    
   const char* dtcuts="sel_BsToJpsif2p_2016MC_UL_Official_Presel_mc.lite_cutopt_bdtr_finalselectioncuts.root";
   const char* dt="sel_BsToJpsif2p_2016MC_UL_Official_Presel_mc.lite_cutopt_bdtr.root";
    

    
/////////////////////////////////////////////////////////////// cuts_101    
    
  //  TCut trig="((JpsiTriggers==1) && dr0<0.1 && dr1<0.1 && mtrkqual==1 && ptrkqual==1)";
  //  TCut resrej1="Q2<11||Q2>8";
  //  TCut f2p="(Phimass>1.525-2.5*0.079 && Phimass<1.525+2.5*0.079)";
  //  TCut b1="(Bmass>5.366-2.5*0.03505 && Bmass<5.366+2.5*0.03505)";
 //   TCut b2="(Bmass>5.6 && Bmass<5.6+5*0.03505 )";
  //  TCut cutm= trig && resrej1 && f2p && b1;
//    TCut cutd= trig && resrej1 && f2p && b2;

/////////////////////////////////////////////////////////////////////////// Followings are the cuts for training purpose   

  //  TCut raretrig="((JpsiTriggers==1 || PsiPTriggers==1 || LMNTTriggers==1) && dr0<0.1 && dr1<0.1 && mtrkqual==1 && ptrkqual==1)";
  //  TCut trig="((JpsiTriggers==1) && dr0<0.1 && dr1<0.1 && mtrkqual==1 && ptrkqual==1)";
  //  TCut resrej1="Q2<11||Q2>8";
  //  TCut resrej2="Q2<8 || Q2>11";
  //  TCut resrej3="Q2<12.5 || Q2>15.0";
  //  TCut f2p="(Phimass>1.525-2.5*0.079 && Phimass<1.525+2.5*0.079)";
  //  TCut b1="(Bmass>5.366-2.5*0.03505 && Bmass<5.366+2.5*0.03505)";
  //  TCut b2="(Bmass>5.6 && Bmass<5.6+5*0.03505)";
  //  TCut cutm= trig && resrej1 && f2p && b1;
  //  TCut cutd= raretrig && f2p && b2 && resrej2 && resrej3;
  
///////////////////////////////////////////// Following  are the final selection cuts


     TCut raretrig="((JpsiTriggers==1 || PsiPTriggers==1 || LMNTTriggers==1) && dr0<0.1 && dr1<0.1 && mtrkqual==1 && ptrkqual==1)";
     TCut resrej1="Q2<11||Q2>8";
     TCut f2p="(Phimass>1.525-2.5*0.079 && Phimass<1.525+2.5*0.079)";
     TCut bdtt="BDT>0.01";
     TCut finalcut = raretrig && resrej1 && f2p && bdtt; 

 

    
    selectSIGNAL_LMNR(dtcuts, path, dt, finalcut);
    
    


}

