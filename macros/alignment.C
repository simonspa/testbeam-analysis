#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include "tools.C"
#include "plotter.C"

using namespace std;

bool cmslogo = false;

void alignment() {
  std::cout << "Run alignment(histogram dir)" << std::endl;
}

void alignment(const char* inputdir, int chip, int startrun, int stoprun) {

  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.2);

  gStyle->SetOptStat(0);

  TGraph *ke = new TGraph(); ke->SetLineColor(kRed+1); ke->SetLineWidth(2);

  TGraph *dx = new TGraph(); dx->SetLineColor(kRed+1); dx->SetLineWidth(3);
  TGraph *dy = new TGraph(); dy->SetLineColor(kBlack); dy->SetLineWidth(3); dy->SetLineStyle(7);
  TGraph *dym = new TGraph(); dym->SetLineColor(kBlack); dym->SetLineWidth(3); dym->SetLineStyle(7);
  TGraph *dz = new TGraph(); dz->SetLineColor(kOrange); dz->SetLineWidth(3);
  TGraph *dtilt = new TGraph(); dtilt->SetLineColor(kRed+1); dtilt->SetLineWidth(2);
  TGraph *dturn = new TGraph(); dturn->SetLineColor(kOrange); dturn->SetLineWidth(2);
  TGraph *drot = new TGraph(); drot->SetLineColor(kOrange); drot->SetLineWidth(2);

  int nruns, nalign = 0, alignrun;
  std::vector<int> telalign;

  bool firstrun = false;
  double nke = 1, ndx = 1, ndy = 1, ndz = 10, ndtilt = 1, ndturn = 1, ndrot = 1;

  double zoffset = 0;
  double xoffset = 0;
  if(chip == 504) zoffset = 50;
  else zoffset = 70;

  if(chip == 506) xoffset = -3;

  // Get all runs for given chip:
  std::vector<int> runs = getruns(inputdir,chip);
  for(std::vector<int>::iterator run = runs.begin(); run != runs.end(); run++) {
    if(*run < startrun || *run > stoprun) continue;

    Int_t alignment = getalignmentrun(inputdir,*run,chip);
    std::vector<Double_t> alignconstants = getalignment(inputdir,*run,chip);

    if(!alignconstants.empty()) {
      // Normalization:
      if(firstrun) {
	firstrun = false;
	nke = alignconstants.at(0);
	ndx = alignconstants.at(1);
	ndy = alignconstants.at(2);
	ndz = alignconstants.at(3);
	ndtilt = alignconstants.at(4);
	ndturn = alignconstants.at(5);
	ndrot = alignconstants.at(6);
      }

      ke->SetPoint(ke->GetN(),ke->GetN(),alignconstants.at(0)/nke);
      dx->SetPoint(dx->GetN(),dx->GetN(),(alignconstants.at(1)-xoffset)/ndx);
      dy->SetPoint(dy->GetN(),dy->GetN(),alignconstants.at(2)/ndy);
      dym->SetPoint(dy->GetN(),dy->GetN(),-1*alignconstants.at(2)/ndy);
      dz->SetPoint(dz->GetN(),dz->GetN(),(alignconstants.at(3)-zoffset)/ndz);
      dtilt->SetPoint(dtilt->GetN(),dtilt->GetN(),alignconstants.at(4)/ndtilt);
      dturn->SetPoint(dturn->GetN(),dturn->GetN(),alignconstants.at(5)/ndturn);
      drot->SetPoint(drot->GetN(),drot->GetN(),alignconstants.at(6)/ndrot);
    }
    if(alignment != alignrun) { nalign++; alignrun = alignment; telalign.push_back(ke->GetN()); }
    cout << "run" << *run << " (chip" << chip << ") align " << alignment << endl;
    nruns++;
  }

  cout << nruns << " runs analyzed with " << nalign << " alignment runs." << endl;

  TString fileName;
  fileName.Form("chip%i-alignment.root",chip);
  TFile * out = TFile::Open(fileName,"RECREATE");
  gDirectory->pwd();

  TCanvas *c1 = new TCanvas();
  TLegend *leg = new TLegend();
  setLegendStyle(leg);
  leg->SetFillColor(kWhite);

  c1->cd();
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";run; alignment shift #left[#mum#right]");
  mg->Add(dx,"l");
  leg->AddEntry(dx, "DUT dx",  "l");

  mg->Add(dy,"l");
  leg->AddEntry(dy, "DUT dy",  "l");

  mg->Add(dz,"l");
  leg->AddEntry(dz, "DUT dz / 10",  "l");

  //mg->Add(dturn,"lp");
  //leg->AddEntry(dturn, "turn",  "lp");

  //mg->Add(drot,"lp");
  //leg->AddEntry(drot,"rotation",  "lp");
  mg->Draw("a");
  if(chip == 506) mg->GetYaxis()->SetRangeUser(-5,5);
  else mg->GetYaxis()->SetRangeUser(-2,2);
  mg->GetXaxis()->SetLimits(-2,dy->GetN()+1);
  DrawFreeCMSLabels(Form("CMS SCM %i  (%1.1f GeV)",chip,5.6),0.045);
  leg->Draw();

  for(int i = 0; i < telalign.size(); i++) {
    TLine *l;
    if(chip == 506) l = new TLine(telalign.at(i)-1,-5,telalign.at(i)-1,5);
    else l = new TLine(telalign.at(i)-1,-2,telalign.at(i)-1,2);
    l->SetLineWidth(2);
    l->SetLineStyle(3);
    l->SetLineColor(kGray+1);
    l->Draw();
  }

  c1->Write();

  // revert some settings:                                                                                                                                   
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetTitleYOffset(1.2);
  gStyle->SetTitleXOffset(1.2);


  TCanvas *c2 = new TCanvas();
  c2->cd();
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad2->SetFrameFillStyle(0);
  pad1->Draw();
  pad1->cd();

  TLegend *leg2 = new TLegend();
  setLegendStyle(leg2);
  leg2->SetFillColor(kWhite);
  dym->Draw("al");
  if(chip == 506) dym->GetYaxis()->SetRangeUser(-3,3);
  else dym->GetYaxis()->SetRangeUser(-1,1);
  dym->GetYaxis()->SetTitle("DUT dy #left[#mum#right]");
  dym->GetXaxis()->SetLimits(-2,dym->GetN()+1);
  leg2->AddEntry(dym, "DUT dy (inverted)",  "l");

  TH1F *h2 = dtilt->GetHistogram();
  TAxis *Ay2 = h2->GetYaxis();
  int maxangle;
  if(chip == 504) maxangle = 60;
  else maxangle = 90;

  Ay2->SetRangeUser(-maxangle, maxangle);
  h2->SetYTitle("DUT #alpha [#circ]");
  pad2->Draw();
  pad2->cd();
  dtilt->GetXaxis()->SetLimits(-2,dym->GetN()+1);
  dtilt->Draw("LAY+");

  leg2->AddEntry(dtilt, "DUT #alpha",  "l");

  DrawFreeCMSLabels(Form("CMS SCM %i  (%1.1f GeV)",chip,5.6),0.045);
  leg2->Draw();

  for(int i = 0; i < telalign.size(); i++) {
    TLine *l = new TLine(telalign.at(i)-1,-maxangle,telalign.at(i)-1,maxangle);
    l->SetLineWidth(2);
    l->SetLineStyle(3);
    l->SetLineColor(kGray+1);
    l->Draw();
  }

  c2->Write();
}
