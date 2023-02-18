#include <TLatex.h>
#include <TString.h>
#include <TLegend.h>
#include <RooCBShape.h>
#include <string>

void cbshapef2p()
{
	gROOT->SetBatch(1);
	gSystem->Load("libRooFit");
	using namespace RooFit;

	const char* path = "/eos/user/t/tthakare/projectWork/exe_files";
	const char* filename = "sel_BsToJpsif2p_2016MC_UL_Official_Presel_mc.lite_cutopt_bdtr_finalselectioncuts.root";

	TChain *tr = new TChain("tree");
	tr->Add(Form("%s/%s", path, filename));

	int nentries_ = tr->GetEntries();
	cout << "\n=> total entries in signal tree = " << nentries_ << endl;

	double bm_min(5.27), bm_max(5.6);
	//double phi_min(1.01), phi_max(1.03);

	RooRealVar Bmass("Bmass","#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}", bm_min, bm_max);
	//RooRealVar Q2("Q2","Q2", 0., 20.);
	//RooRealVar Phimass("Phimass", "#bf{#phi}", phi_min, phi_max);
	RooArgSet  observables(Bmass);

	//create dataset from the input file
	RooDataSet data("data","dataset with Bmass", tr, observables);

	TCut c1 = Form("Bmass>%f && Bmass<%f",bm_min,bm_max);
	//TCut ph = Form("Phimass>%f && Phimass<%f", phi_min, phi_max);
        TCut cutTotal = c1; //&& ph;

	//reduced data set- apply the Bmass cut(fitting range)
	RooDataSet *redData = (RooDataSet*)data.reduce(cutTotal);
	std::cout<<"After final cut: "<<redData->sumEntries()<<std::endl;

	//define double crystal ball pdf for fitting
	RooRealVar  mean("mean","common means for Crystal Balls", 5.367, bm_min, bm_max);
	RooRealVar  sigma1("sigma1","sigma of CB1",  0.024, 0.002, 0.031);
     RooRealVar  sigma2("sigma2","sigma of CB2",  0.045, 0.001, 0.1);
     RooRealVar  sigM_frac("sigM_frac","fraction of CB", 0.44, 0., 1.);
     RooRealVar  n1("n1", "", 0.5);
     RooRealVar  n2("n2", "", 10.);
     RooRealVar  alpha1("alpha1","alpha for CB1", 1., 1., 2.);
     RooRealVar  alpha2("alpha2","alpha for CB2", -1., -2.0, -1.);

	RooCBShape CB1("CB1","Crystal Ball-1", Bmass, mean, sigma1, alpha1,n1);
        RooCBShape CB2("CB2","Crystal Ball-2", Bmass, mean, sigma2, alpha2,n2);

	RooAddPdf CB("CB","CB1+CB2", RooArgList(CB1,CB2), RooArgList(sigM_frac));
        RooRealVar nsig("nsig","nsig",1E3,500,1.1*tr->GetEntries());

	 //final model used for fitting
        RooExtendPdf model("model","model", CB, nsig);

	model.Print();

	cout<< "pdf evaluation for MC fitting" << endl;
        cout<< model.getLogVal() <<endl;

	RooFitResult* fitres = model.fitTo(*redData, Minos(true), Extended(true), Save(true));

	TCanvas *c = new TCanvas("c","c",800, 700);
        TPad *p1   = new TPad("p1","p1", 0.01, 0.25, 0.995, 0.97);
        TPad *p2   = new TPad("p2","p2", 0.01, 0.02, 0.995, 0.24);
        p1->Draw();
        p2->Draw();

        p1->cd();
  //      gPad->SetLogy();

	// creating frame for Bmass invariant mass distribution
        RooPlot *xframe = Bmass.frame(Title(""), Bins(400));
        xframe->SetTitle("");
        xframe->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");

	// plotting data of the Bmass from file on the frame
        redData->plotOn(xframe,RooFit::Name("data"));

	// plotting final pdf on the frame
        model.plotOn(xframe, RooFit::Name("Full PDF"), LineColor(kBlue));

	// create a pull distribution for the fit
        RooHist* hpull = xframe->pullHist();

	printf("**************************************");
        printf("\n #sigma_{1}: %.7f", sigma1.getVal());
        printf("\n #sigma_{2}: %.7f", sigma2.getVal());
        double eff_sigma = sqrt(sigM_frac.getVal()*pow(sigma1.getVal(),2) + (1 - sigM_frac.getVal()) * pow(sigma2.getVal(),2));
        printf("\n #sigma_{eff}: %.7f", eff_sigma);
        printf("\n************************************\n");

	int nFloatParam = fitres->floatParsFinal().getSize();
	cout << "number of floating parameters-> " << nFloatParam << endl;
	double chi2dof = xframe->chiSquare(nFloatParam);
        std::cout<<"\n"<<std::endl;
        std::cout<<"#chi^{2}/dof= "<< chi2dof << std::endl;

	// plotting individual crystal ball in the frame
        model.plotOn(xframe, RooFit::Name("CB1"), Components("CB1"), LineStyle(kDashed), LineColor(kRed));
        model.plotOn(xframe, RooFit::Name("CB2"), Components("CB2"), LineStyle(kDashed), LineColor(kGreen));

	xframe->Draw();

	TPaveText* paveText = new TPaveText(0.65,0.40,0.83,0.88,"NDC");
	paveText->SetBorderSize(0.0);
        paveText->SetFillColor(kWhite);
        paveText->SetFillStyle(0);
        paveText->SetTextSize(0.02);
	paveText->AddText(Form("N_{sig} = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));
	paveText->AddText(Form("#mu  = %.7f #pm %.7f GeV" , mean.getVal() , mean.getError()));
	//paveText->AddText(Form("#sigma_{1} = %.7f #pm %.7f GeV", sigma1.getVal(), sigma1.getError()));
	//paveText->AddText(Form("#sigma_{2} = %.7f #pm %.7f GeV", sigma2.getVal(), sigma2.getError()));
	paveText->AddText(Form("#sigma_{eff} = %.7f  GeV", eff_sigma));
	paveText->AddText(Form("f = %.7f #pm %.7f ", sigM_frac.getVal(), sigM_frac.getError()));
	paveText->AddText(Form("#alpha_{1} = %.7f #pm %.7f ", alpha1.getVal(), alpha1.getError()));
	paveText->AddText(Form("#alpha_{2} = %.7f #pm %.7f ", alpha2.getVal(), alpha2.getError()));
	//paveText->AddText(Form("n_{1} = %.7f #pm %.7f", n1.getVal(), n1.getError()));
        //paveText->AddText(Form("n_{2} = %.7f #pm %.7f", n2.getVal(), n2.getError()));
        paveText->AddText(Form("#chi^{2}/dof  = %.5f ", chi2dof ));

	/*if (fitres != NULL){
    if ((fitres->covQual() == 3) && (fitres->status() == 0)){
      if (paveText != NULL) paveText->AddText("Fit Status: GOOD");
    } else {
      if (paveText != NULL) paveText->AddText("Fit Status: BAD");
    }
  }*/

	paveText->Draw();

	TLatex *mark = new TLatex();
        mark->SetNDC(true);
        double startY = 0.92;
        mark->SetTextFont(42);
        mark->SetTextSize(0.035);
        mark->DrawLatex(0.12,startY,"#bf{CMS} #it{Simulation}");
        mark->DrawLatex(0.83,startY, "#it{13 TeV}");

    //    mark->DrawLatex(0.75,startY,"#scale[0.8]{66226.56 fb^{-1} (13 TeV)}");

        mark->Draw();

	p2->cd();

        RooPlot* pull = Bmass.frame(Title(""));
        pull->addPlotable(hpull, "BX");
        pull->SetMinimum(-5);
        pull->SetMaximum(5);
        hpull->SetFillColor(kBlue);
        pull->GetYaxis()->SetNdivisions(5);
        pull->SetTitle("#bf{Pull Distribution}");
        pull->GetXaxis()->SetTitle("#bf{m(K^{+}K^{-}#mu^{+}#mu^{-}) [GeV]}");
        pull->GetYaxis()->SetTitle("");
        pull->GetYaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull->GetXaxis()->SetLabelSize( (0.7/0.3)*xframe->GetXaxis()->GetLabelSize() );
        pull->GetXaxis()->SetTitle("");
        pull->GetXaxis()->SetLimits(bm_min, bm_max);
        pull->Draw("P");

        TLine* l1 = new TLine(5.27, 3, 5.6, 3);
        TLine* l2 = new TLine(5.27, -3, 5.6, -3);
        l1->SetLineColor(2);
        l2->SetLineColor(2);
        TLine* l3 = new TLine(5.27, 0, 5.6, 0);
        l3->SetLineColor(2);
        l1->Draw();
        l2->Draw();
        l3->Draw();

	c->SaveAs(Form("%s/cb_final_sel.pdf", path));
	c->SaveAs(Form("%s/cb_final_sel.png", path));

}







