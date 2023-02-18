#include<iostream>
using namespace std;

void plot(TH1D *h1, TH1D *h2, int i,const char* varname, const char* title)
{
THStack *hs = new THStack("hs", "stacked histograms");


h1->SetStats(0);
h2->SetStats(0);

h1->SetLineWidth(2);
h2->SetLineWidth(2);

h1->SetLineColor(2);
h2->SetLineColor(4);

h1->SetFillColor(2);
h2->SetFillColor(4);

h1->SetFillStyle(3004);
h2->SetFillStyle(3001);

hs->Add(h1,"hist1");
hs->Add(h2,"hist2");
TCanvas *c1 = new TCanvas("Final_amc", "Application Plots ", 800, 700);

hs->Draw("no stack");
hs->SetTitle("");
hs->GetXaxis()->SetTitle(title);
hs->GetYaxis()->SetTitle("#bf{A.U.}");

TLegend *leg=new TLegend(.1,.8,.2,.9,"Lables");
leg->SetFillColor(0);
leg->AddEntry(h1,"Before BDT");
leg->AddEntry(h2,"After BDT");
leg->DrawClone("Same");

c1->SaveAs(Form("plot_final_sel_cut_%i_%s.png", i,varname));
delete c1;
}
int stackplotploter()
{
const int len =5;
float xmin[len]={5.27,-0.4, 0, 0,-1};
float xmax[len]={5.8,0.4, 70, 1, 1};
const char* varname[len]={"Bmass ","BDT","Blxysig", "Bvtxcl","BDTG"};
const char* title[len]={"m(K^{+}K^{-}#mu^{+}#mu^{-}) GeV","BDT","B_{s}L_{xy}/#sigma","B_{s} vtx. CL","BDTG"};

TH1D *hi1;
TH1D *hi2;

TChain *ch1 = new TChain("tree");
ch1->Add("sel_Combine_2016_Mini_Presel_data_cutopt_bdtr.root");
TTree *tr1 = ch1;

TChain *ch2 = new TChain("tree");
ch2->Add("sel_Combine_2016_Mini_Presel_data_cutopt_bdtr_finalselectioncuts.root");
TTree *tr2 = ch2;

int norm =1;
double sc1, sc2;

for(int i=0;i<len;i++)
{

hi1=new TH1D("hi1","hist1",80,xmin[i],xmax[i]);
hi2=new TH1D("hi2","hist2",80,xmin[i],xmax[i]);

tr1->Project("hi1", varname[i]);
tr2->Project("hi2", varname[i]);

sc1= norm/hi1->Integral();
sc2= norm/hi2->Integral();

 hi1->Scale(sc1);
 hi2->Scale(sc2);

plot(hi1,hi2, i, varname[i],title[i]);

delete hi1;
delete hi2;
}

return 0;
}
