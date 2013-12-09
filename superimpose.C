#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TPaveStats.h"
#include "TStyle.h"

#include <iostream>
#include <string>

using namespace std;

TString getRightAxis(TString);
double whichRangeLow (TString, TString);
double whichRangeHigh (TString, TString);
int howRebin(TString);
TLegend*  MakeLegend(void);

void superimpose(TString HISTO, TString SCENARIO, TString REGION, TString DIRS) {

  //TString SCENARIO = "0100";  //this fixes sample
  //TString DIRS = "GLB"; //FMS //this fixes dir
  //TString HISTO = "pt";
  //TString REGION = "barrel24d";

  //passing cuts
  TFile *f = new TFile("singleMu"+SCENARIO+".root");
  f->cd();
  gDirectory->cd(DIRS+"segsInFits");
  //get histograms
  //knopwn histonames: phi, pt, p, eta
  //known whats: endcap, barrel, barrel24d, barrel4d, overlap

  cout << "segs_"+REGION+"_pT"+SCENARIO+"_"+HISTO << endl;

  TH1D* hseg = (TH1D*)gDirectory->Get("segs_"+REGION+"_pT"+SCENARIO+"_"+HISTO);
  TH1D* hhit = (TH1D*)gDirectory->Get("hits_"+REGION+"_pT"+SCENARIO+"_"+HISTO);

   gStyle->SetEndErrorSize(2);
   gStyle->SetErrorX(0.5);

  //draw
  TCanvas* c = new TCanvas();
  c->cd();
  //c->SetLogx();
  hseg->GetXaxis()->SetTitle(getRightAxis(HISTO));
  hseg->GetYaxis()->SetTitle("Entries");
  hseg->SetLineColor(kRed);
  hseg->SetLineWidth(1.2);

  if (HISTO == "pt" || HISTO == "p") hseg->Rebin(howRebin(SCENARIO));
  hseg->GetXaxis()->SetRangeUser(whichRangeLow(SCENARIO,REGION),whichRangeHigh(SCENARIO,REGION));
  hseg->Draw("E");  

  hhit->SetLineWidth(1.2);
  if (HISTO == "pt" || HISTO == "p") hhit->Rebin(howRebin(SCENARIO));
  hhit->Draw("Esames"); 
  c->Update();

  TPaveStats *p1 = (TPaveStats*)hseg->GetListOfFunctions()->FindObject("stats");
  hseg->GetListOfFunctions()->Remove(p1);
  hseg->SetStats(0);
  p1->SetX1NDC(.66);
  p1->SetX2NDC(.82);
  p1->SetTextColor(kRed);
  p1->Draw();

  TLegend *leg = MakeLegend();
  leg->AddEntry(hseg,"Segment fit","l");
  leg->AddEntry(hhit,"RecHit fit","l");
  leg->Draw("same");
  
  c->SaveAs(DIRS+SCENARIO+REGION+"_"+HISTO+".png");

}

TString getRightAxis(TString what) {

  if (what == "pt") return "p_{T} [GeV]";
  else if (what == "p") return "p [GeV]";
  else if (what == "eta") return "#eta";
  else if (what == "phi") return "#phi";
  else return "";
}

double whichRangeLow (TString sample, TString region) {

  if (region == "barrel" || region == "overlap") {
     if (sample == "1000") return 500;
     else if (sample == "0100") return 85;
     else if (sample == "0010") return 8;
  } else if (region == "endcap") {
     if (sample == "1000") return 200;
     else if (sample == "0100") return 70;
     else if (sample == "0010") return 7;
  }

  return 0;
}

double whichRangeHigh (TString sample, TString region) {

  if (region == "barrel" || region == "overlap") {
     if (sample == "1000") return 1500;
     else if (sample == "0100") return 115;
     else if (sample == "0010") return 12;
  } else if (region == "endcap") {
     if (sample == "1000") return 1800;
     else if (sample == "0100") return 130;
     else if (sample == "0010") return 13;
  }

  return 0;
}


int howRebin(TString sample) {

  if (sample == "1000") return 15;
  else if (sample == "0100") return 5;
  else if (sample == "0010") return 2;
  else return 1;

}

TLegend*  MakeLegend() {

   TLegend *_leg= new TLegend(0.75,0.45,0.90,0.62);
   _leg->SetTextFont(72);
   _leg->SetTextSize(0.035);
   _leg->SetBorderSize(0);
   _leg->SetFillColor(0);
   return _leg;

}

