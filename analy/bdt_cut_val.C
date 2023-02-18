#include <sstream>

void bdt_cut_val(){
    
    
  const std::string filenameS("./sel_BsToJpsif2p_2016MC_UL_Official_Presel_mc.lite_cutopt_bdtr.root");
  const std::string trnameS = "tree";
  const std::string filenameB("./sel_Combine_2016_Mini_Presel_data_cutopt_bdtr.root");
  const std::string trnameB = "tree";
    

    
  TFile* fileS = TFile::Open( filenameS.c_str() );
  if( !fileS ) std::cout << "file " << filenameS << " does not exist" << std::endl;
  TTree* tree = (TTree*)fileS->Get( trnameS.c_str() );
  if( !tree ) std::cout << "tree " << trnameS << " does not exist" << std::endl;
    

  TFile* fileB = TFile::Open( filenameB.c_str() );
  if( !fileB ) std::cout << "file " << filenameB << " does not exist" << std::endl;
  TTree* tree1 = (TTree*)fileB->Get( trnameB.c_str() );
  if( !tree1 ) std::cout << "tree " << trnameB << " does not exist" << std::endl;
    
    
  const Int_t nBins = 200 ; 
    
  double bdt_cuts[nBins];
  double efficiencies2[nBins];
  double efficiencies2_error[nBins];
 
        
  double effvals[nBins];
    
  //********************
  //loop starting here
  //********************
  for(int i=0; i < nBins; i=i+1) {
    double step_size = double(2./nBins);
    cout << "step size = " << step_size << endl;
   

    double cut_val = -1.0 + step_size*i ;
    
    bdt_cuts[i] = cut_val;

    std::cout << i << std::endl;        
    std::stringstream c;

      c << "BDT" << " >= " << cut_val<<"&& ((JpsiTriggers==1) && dr0<0.1 && dr1<0.1 && mtrkqual==1 && ptrkqual==1) && (Q2<11||Q2>8) && ((Phimass>1.525-2.5*0.079 && Phimass<1.525+2.5*0.079)) && ((Bmass>5.366-2.5*0.03505 && Bmass<5.366+2.5*0.03505))";
    const std::string cut = c.str();
        
    std::cout << cut << std::endl;

    // additional set of cuts to pick the background                                                                                      
    std::stringstream r;
    
    r << "BDT" << " >= " << cut_val << " && ((Bmass>5.6 && Bmass<5.6+5*0.03505)) && ((JpsiTriggers==1 || PsiPTriggers==1 || LMNTTriggers==1)) && (dr0<0.1) && (dr1<0.1) && (mtrkqual==1) &&( ptrkqual==1) && ((Phimass>1.525-2.5*0.079 && Phimass<1.525+2.5*0.079)) && (Q2<8 || Q2>11) && (Q2<12.5 || Q2>15.0)";
    const std::string cut2 = r.str();

    std::cout << cut2 << std::endl;

  
   
    double integB_pre = tree1->GetEntries();
    
    double integB     = tree1->GetEntries(cut2.c_str());
    

    double MC_pre  = tree->GetEntries();
    double MC_post = tree->GetEntries(cut.c_str());

    double SF = 35.9/2788.4;
    double MC_post_scaled = MC_post*SF ;

    double efficiency2 = (MC_post_scaled)/(sqrt(MC_post_scaled + integB));
    if (MC_post_scaled == 0 && integB == 0) efficiency2 = 0; 
    efficiencies2[i] = efficiency2;

    std::cout << "B = " << integB << std::endl;
    std::cout << "S = " << MC_post << std::endl;
    std::cout << "Scale Factor = " << SF << std::endl;
    std::cout << "S_scaled = " << MC_post << "*" << SF << " = " << MC_post_scaled << std::endl;
    std::cout << "Significance = S/sqrt(S+B) = " << efficiency2 << std::endl;
    std::cout << "++++++++++++++++++++++++++++" << std::endl;

  }

    //The maximum significance value 
    double *i1;
    i1 = std::max_element(efficiencies2, efficiencies2 + nBins); 
    std::cout << "max significance value  = " << *i1 << std::endl;


    double max_sign = -99;
    int maxCut = 99;
    for (int m = 0 ; m < nBins ; m++){
      if (efficiencies2[m] > max_sign){
	max_sign = efficiencies2[m];
	maxCut = m;
      }
    }

    std::cout << "max. significance value = " << max_sign << " for BDT cut value = " << bdt_cuts[maxCut] << std::endl;


    TCanvas *c2 = new TCanvas("c2", "",800,600);
    TGraph* graph2 = new TGraph(nBins, bdt_cuts, efficiencies2);
    graph2->SetTitle("FOM");
    graph2->GetYaxis()->SetLabelSize(0.05);
    graph2->SetMarkerStyle(20);
    graph2->GetXaxis()->SetTitle("BDT cut > ");  
    graph2->GetXaxis()->SetRangeUser(-1.0,1.0);
    graph2->GetYaxis()->SetTitle("S/sqrt(S+B)");
    
    graph2->Draw("APL");
    c2->SaveAs("fom.png");    
}








