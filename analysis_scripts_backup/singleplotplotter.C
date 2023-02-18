#include<iostream>
using namespace std;

void plot(TH1D *h1, int i,const char* varname, const char* title)
{



h1->SetStats(0);


h1->SetLineWidth(2);

h1->SetLineColor(4);


h1->SetFillColor(4);


h1->SetFillStyle(3001);



TCanvas *c1 = new TCanvas("Final_amc", "Application Plots ", 800, 700);

h1->Draw();
h1->SetTitle("");
h1->GetXaxis()->SetTitle(title);
h1->GetYaxis()->SetTitle("#bf{A.U.}");

TLegend *leg=new TLegend(.1,.8,.2,.9,"Lables");
leg->SetFillColor(0);
leg->AddEntry(h1,"After BDT");

leg->DrawClone("Same");

c1->SaveAs(Form("single_plot_%i_%s.png", i,varname));
delete c1;
}
int singleplotplotter()
{
const int len =5;
float xmin[len]={5.27,-0.4, 0, 0,-1};
float xmax[len]={5.8,0.4, 70, 1, 1};
const char* varname[len]={"Bmass ","BDT","Blxysig", "Bvtxcl","BDTG"};
const char* title[len]={"m(K^{+}K^{-}#mu^{+}#mu^{-}) GeV","BDT","B_{s}L_{xy}/#sigma","B_{s} vtx. CL","BDTG"};

TH1D *hi1;


TChain *ch1 = new TChain("tree");
ch1->Add("sel_Combine_2016_Mini_Presel_data_cutopt_bdtr_finalselectioncuts.root");
TTree *tr1 = ch1;

TChain *ch2 = new TChain("tree");
ch2->Add("sel_Combine_2016_Mini_Presel_data_cutopt_bdtr_finalselectioncuts.root");
TTree *tr2 = ch2;

int norm =1;
double sc1, sc2;

for(int i=0;i<len;i++)
{

hi1=new TH1D("hi1","hist1",80,xmin[i],xmax[i]);


tr1->Project("hi1", varname[i]);

//sc1= norm/hi1->Integral();


 //hi1->Scale(sc1);
 

plot(hi1, i, varname[i],title[i]);

delete hi1;

}

return 0;
}
