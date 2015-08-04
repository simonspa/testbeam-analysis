#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include "tools.C"

using namespace std;

void nrows() {
  std::cout << "Run nrows(histogram dir)" << std::endl;
}

void nrows(const char* inputdir, int chip, int startrun, int stoprun) {

  TCanvas *c1 = new TCanvas("c1","nrows",600,600);
  TProfile *nrows = new TProfile("nrows"," ",130,5,85,0,60,"");

  TCanvas *c2 = new TCanvas("c2","nrowstan",600,600);
  TProfile *nrowstan = new TProfile("nrowstan"," ",130,0,12,0,60,"");

  gStyle->SetOptStat(0);

  int nruns, nclusters;

  // Get all runs for given chip:
  std::vector<int> runs = getruns(inputdir,chip);
  for(std::vector<int>::iterator run = runs.begin(); run != runs.end(); run++) {
    if(*run < startrun || *run > stoprun) continue;

    TString fileName;
    fileName += inputdir;
    if( !fileName.EndsWith("/") ) fileName += "/";
    fileName += "run" + ZeroPadNumber(*run,6) + "-analysis.root";

    TFile *source;
    if (!gSystem->AccessPathName( fileName )) source = TFile::Open(fileName);
    else continue;

    //cout << "File: " << fileName << endl;

    source->cd("MyEUTelAnalysisCMSPixel");
    TH1 *h;
    gDirectory->GetObject("cmsnrow",h);
    if(!h) continue;

    Double_t nrow = h->GetMean();
    Double_t tilt = gettilt(inputdir,*run,chip);

    // Collect statistics:
    nclusters += ((TH1D*)gDirectory->Get("cmsnrow"))->GetEntries();

    cout << "run" << *run << " (chip" << chip << ") nrow " << nrow << " tilt " << tilt << endl;
    nrows->Fill(tilt,nrow,1);
    nrowstan->Fill(TMath::Tan(tilt/180*TMath::Pi()),nrow,1);

    nruns++;
    delete source;
  }

  cout << nruns << " runs analyzed with " << nclusters << " clusters in total." << endl;

  TLegend * leg = new TLegend();

  c1->cd();
  nrows->SetTitle(";tilt angle [#deg];rows per cluster");
  nrows->SetMarkerStyle(20);
  nrows->SetMarkerColor(1);
  nrows->Draw();

  std::vector<double> vtilt = getsimulation("tilt", chip);
  std::vector<double> vtilttan = getsimulation("tilttan", chip);
  std::vector<double> vnrow = getsimulation("ncol", chip);

  TGraph *si = new TGraph( vtilt.size(), &(vtilt[0]), &(vnrow[0]) ); // sim
  si->SetLineColor(2);
  si->SetLineWidth(3);
  si->SetMarkerColor(2);
  nrows->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
  si->Draw("PL"); // without axis option: overlay

  leg->Draw();

  c2->cd();
  nrowstan->SetTitle(";tan(#alpha);rows per cluster");
  nrowstan->SetMarkerStyle(20);
  nrowstan->SetMarkerColor(1);
  nrowstan->Draw();

  TGraph *si_tan = new TGraph( vtilttan.size(), &(vtilttan[0]), &(vnrow[0]) ); // sim
  si_tan->SetLineColor(2);
  si_tan->SetLineWidth(3);
  si_tan->SetMarkerColor(2);
  nrowstan->GetXaxis()->SetRangeUser(vtilttan.front(), vtilttan.back());
  si_tan->Draw("PL"); // without axis option: overlay

}
