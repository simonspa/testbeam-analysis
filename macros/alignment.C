#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
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
  gStyle->SetTitleYOffset(1);

  gStyle->SetOptStat(0);

  TGraph *ke = new TGraph(); ke->SetLineColor(kRed+1);

  TGraph *dx = new TGraph(); dx->SetLineColor(kBlue);
  TGraph *dy = new TGraph(); dy->SetLineColor(kRed);
  TGraph *dz = new TGraph(); dz->SetLineColor(kGreen);
  TGraph *dtilt = new TGraph(); dtilt->SetLineColor(kOrange);
  TGraph *dturn = new TGraph(); dturn->SetLineColor(kOrange);
  TGraph *drot = new TGraph(); drot->SetLineColor(kOrange);

  int nruns, nalign = 0, alignrun;

  // Get all runs for given chip:
  std::vector<int> runs = getruns(inputdir,chip);
  for(std::vector<int>::iterator run = runs.begin(); run != runs.end(); run++) {
    if(*run < startrun || *run > stoprun) continue;

    Int_t alignment = getalignmentrun(inputdir,*run,chip);
    std::vector<Double_t> alignconstants = getalignment(inputdir,*run,chip);

    if(!alignconstants.empty()) {
      ke->SetPoint(ke->GetN(),ke->GetN(),alignconstants.at(0));
      dx->SetPoint(dx->GetN(),dx->GetN(),alignconstants.at(1));
      dy->SetPoint(dy->GetN(),dy->GetN(),alignconstants.at(2));
      dz->SetPoint(dz->GetN(),dz->GetN(),alignconstants.at(3));
      dtilt->SetPoint(dtilt->GetN(),dtilt->GetN(),alignconstants.at(4));
      dturn->SetPoint(dturn->GetN(),dturn->GetN(),alignconstants.at(5));
      drot->SetPoint(drot->GetN(),drot->GetN(),alignconstants.at(6));
    }
    if(alignment != alignrun) { nalign++; alignrun = alignment;}
    cout << "run" << *run << " (chip" << chip << ") align " << alignment << endl;
    nruns++;
  }

  cout << nruns << " runs analyzed with " << nalign << " alignment runs." << endl;

  TCanvas *c1 = new TCanvas();
  c1->cd();
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(dx,"lp");
  mg->Add(dy,"lp");
  //mg->Add(dz,"lp");
  mg->Add(dtilt,"lp");
  mg->Add(dturn,"lp");
  mg->Add(drot,"lp");
  mg->Draw("a");

  TCanvas *c2 = new TCanvas();
  c2->cd();
  ke->Draw("apl");
}
