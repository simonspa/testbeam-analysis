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

Double_t restel = 4.8;

void ncols() {
  std::cout << "Run ncols(histogram dir)" << std::endl;
}

void ncols(const char* inputdir) {
  ncols(inputdir,506,0,99999);
}

void ncols(const char* inputdir, int chip, int startrun, int stoprun) {

  TCanvas *c1 = new TCanvas("c1","ncols",600,600);
  TProfile *ncols = new TProfile("ncols"," ",130,5,85,0,60,"");

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
    gDirectory->GetObject("cmsncol",h);
    if(!h) continue;

    Double_t ncol = h->GetMean();
    Double_t tilt = gettilt(inputdir,*run,chip);

    // Collect statistics:
    nclusters += ((TH1D*)gDirectory->Get("cmsncol"))->GetEntries();

    cout << "run" << *run << " (chip" << chip << ") ncol " << ncol << " tilt " << tilt << endl;
    ncols->Fill(tilt,ncol,1);

    nruns++;
    delete source;
  }

  cout << nruns << " runs analyzed with " << nclusters << " clusters in total." << endl;

  TLegend * leg = new TLegend();

  c1->cd();
  ncols->SetTitle(";tilt angle [#deg];columns per cluster");
  ncols->SetMarkerStyle(20);
  ncols->SetMarkerColor(1);
  ncols->Draw();

  std::vector<double> vtilt = getsimulation("tilt", chip);
  std::vector<double> vncol = getsimulation("ncol", chip);

  TGraph *si = new TGraph( vtilt.size(), &(vtilt[0]), &(vncol[0]) ); // sim
  si->SetLineColor(2);
  si->SetLineWidth(3);
  si->SetMarkerColor(2);
  si->Draw("PL"); // without axis option: overlay

  leg->Draw();
}
