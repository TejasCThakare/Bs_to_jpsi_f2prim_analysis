#include<iostream>
using namespace std;

void plot(TH1D *h1, TH1D *h2, int i, const char* varname, const char* title)
{
	gROOT->SetBatch(1);
THStack *hs = new THStack("hs", "stacked histograms");


h1->SetStats(0);
h2->SetStats(0);

h1->SetLineWidth(2);
h2->SetLineWidth(2);

h1->SetLineColor(4);
h2->SetLineColor(2);

h1->SetFillColor(4);
h2->SetFillColor(2);

h1->SetFillStyle(3001);
h2->SetFillStyle(3004);

hs->Add(h1,"hist1");
hs->Add(h2,"hist2");

TCanvas *c1 = new TCanvas("Final_amc", "Analysis Plot", 800, 700);

hs->Draw("no stack");
hs->SetTitle("");
hs->GetXaxis()->SetTitle(title);
hs->GetYaxis()->SetTitle("#bf{A.U.}");

TLegend *leg=new TLegend(.6,.7,.9,.9,"");
leg->SetFillColor(0);
leg->AddEntry(h1,"MC");
leg->AddEntry(h2,"Data");
leg->DrawClone("Same");

c1->SaveAs(Form("plot_%i.png", i));
delete c1;
}
int plotamcf2()
{
//cuts
//float sigma=0.037;
TCut trig="((JpsiTriggers==1) && dr0<0.1 && dr1<0.1 && mtrkqual==1 && ptrkqual==1)";
TCut resrej1="Q2<11||Q2>8";
TCut f2p="(Phimass>1.525-2.5*0.079 && Phimass<1.525+2.5*0.079)";
TCut b1="(Bmass>5.366-2.5*0.03505 && Bmass<5.366+2.5*0.03505)";
TCut b2="(Bmass>5.6)";
TCut cutm= trig && resrej1 && f2p && b1;
TCut cutd= trig && resrej1 && f2p && b2;


const int len =42;
float xmin[len]={0,0,0,0,0,0,0,0,0,0,0.98,0,0,0.6,5.2,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.5,0,0,-3.3,-1.1,0,0,0,0,-1.1,-1.1,-3.3};
float xmax[len]={14,0.01,0.01,30,30,0.01,0.01,14,14,70,1.001,1.001,6,1.001,6.32,55,1.001,1.001,20,20,25,25,35,35,0.02,0.02,0.007,0.007,1.001,1.001,2.5,3.2,20,3.3,1.1,40,1.1,65,0.05,1.1,1.1,3.3};



const char* varname[len]= {"TMath::Max(Kmpt, Kppt)","MumMinIP/MumMinIPE", "MupMinIP/MupMinIPE", "MumMinIPE", "MupMinIPE", "KmtrkMinIP/KmtrkMinIPE", "KptrkMinIP/KptrkMinIPE", "KmtrkMinIPE", "KptrkMinIPE", "Blxysig", "Bcosalphabs2d", "Bvtxcl", "Bsdcasigbs", "BsIso", "Bmass", "Bpt", "kmtrkIso", "kptrkIso", "Kmtrkdcasigbs", "Kptrkdcasigbs", "Mumpt", "Muppt", "Mumdcasigbs", "Mupdcasigbs", "MumMinIP2D/MumMinIP2DE", "MupMinIP2D/MupMinIP2DE", "MumMinIP2DE", "MupMinIP2DE", "mupIso", "mumIso", "Beta", "Bphi", "Phipt", "Phiphi", "Bcosalphabs", "dimupt", "dimuvtxcl", "dimulsig", "dimuDCA", "CosThetaL", "CosThetaK", "Phi"};

const char* title[len] = {"Max(K^{#pm} p_{T}) GeV","#mu^{-} min. IP/#sigma", "#mu^{+} min. IP/#sigma", "#mu^{-} min. IPE", "#mu^{+} min. IPE","k^{-} trk. min. IP/#sigma","k^{+} trk. min. IP/#sigma","K^{-} trk. min. IPE","K^{+} trk. min. IPE" , "B_{s} L_{xy}/#sigma","B_{s} cos#alpha_{xy}","B_{s} vtx. CL","B_{s} DCA/#sigma","B_{s} Iso.","m(K^{+}K^{-}#mu^{+}#mu^{-}) GeV","B_{s} p_{T} GeV","K^{-} trk. Iso.","K^{+} trk. Iso.","K^{-} trk. DCA/#sigma","K^{+} trk. DCA/#sigma","#mu^{-} p_{T} GeV","#mu^{+} p_{T} GeV","#mu^{-} DCA/#sigma","#mu^{+} DCA/#sigma","#mu^{-} min. IP_{xy}/#sigma","#mu^{+} min. IP_{xy}/#sigma","#mu^{-} min. IP2DE","#mu^{+} min. IP2DE","#mu^{+} Iso.","#mu^{-} Iso.","#Beta","B_{s}#phi","#phi p_{T} GeV","#phi#phi","B_{s} cos#alpha ","#mu^{+}#mu^{-} p_{T} GeV","#mu^{+}#mu^{-} vtx. CL","#mu^{+}#mu^{-} L/#sigma","#mu^{+}#mu^{-} DCA.","cos#theta_{L}","cos#theta_{k}","#phi"};

TH1D *hi1;
TH1D *hi2;

TChain *ch1 = new TChain("tree");
ch1->Add("mc.root");
TTree *tr1 = ch1;

TChain *ch2 = new TChain("tree");
ch2->Add("data.root");
TTree *tr2 = ch2;

int norm =1;
double sc1, sc2;

for(int i=0;i<len;i++)
{

hi1=new TH1D("hi1","hist1",80,xmin[i],xmax[i]);
hi2=new TH1D("hi2","hist2",80,xmin[i],xmax[i]);

tr1->Project("hi1", varname[i], cutm);
tr2->Project("hi2", varname[i], cutd);

sc1= norm/hi1->Integral();
sc2= norm/hi2->Integral();

hi1->Scale(sc1);
hi2->Scale(sc2);

plot(hi1,hi2, i,  varname[i],title[i]);

delete hi1;
delete hi2;
}

return 0;
}
