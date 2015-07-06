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

void clustercharge() {
  std::cout << "Run clustercharge(histogram dir)" << std::endl;
}

void clustercharge(const char* inputdir) {
  clustercharge(inputdir,506,0,99999);
}

void clustercharge(const char* inputdir, int chip, int startrun, int stoprun) {

  TCanvas *c1 = new TCanvas("c1","clustercharge",600,600);
  TProfile *clustercharge = new TProfile("clustercharge"," ",130,0,9,0,300,"");

  TCanvas *c2 = new TCanvas("c2","clustercharge",600,600);
  TProfile *clustercharge_tilt = new TProfile("clustercharge_tilt"," ",130,0,85,0,300,"");

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
    gDirectory->GetObject("cmsqf",h);
    if(!h) continue;

    Double_t landau = fitfulllang("cmsqf");
    Double_t tilt = gettilt(inputdir,*run,chip);
    Double_t cos = 1/TMath::Cos(TMath::TwoPi()*tilt/360);

    // Collect statistics:
    nclusters += ((TH1D*)gDirectory->Get("cmsqf"))->GetEntries();

    cout << "run" << *run << " (chip" << chip << ") peak " << landau << " tilt " << tilt << " 1/cos " << cos << endl;
    clustercharge->Fill(cos,landau,1);
    clustercharge_tilt->Fill(tilt,landau,1);

    nruns++;
    delete source;
  }

  cout << nruns << " runs analyzed with " << nclusters << " clusters in total." << endl;

  TLegend * leg = new TLegend();

  c1->cd();
  clustercharge->SetTitle(";1/cos(tilt angle);peak cluster charge");
  clustercharge->SetMarkerStyle(20);
  clustercharge->SetMarkerColor(1);
  clustercharge->Draw();

  std::vector<double> vpath = getsimulation("path", chip);
  std::vector<double> vtilt = getsimulation("tilt", chip);
  std::vector<double> vpeak = getsimulation("peak", chip);

  TGraph *si = new TGraph( vpath.size(), &(vpath[0]), &(vpeak[0]) ); // sim
  si->SetLineColor(2);
  si->SetLineWidth(3);
  si->SetMarkerColor(2);
  si->Draw("PL"); // without axis option: overlay

  leg->Draw();

  c2->cd();
  clustercharge_tilt->SetTitle(";tilt angle [deg];peak cluster charge");
  clustercharge_tilt->SetMarkerStyle(20);
  clustercharge_tilt->SetMarkerColor(1);
  clustercharge_tilt->Draw();

  TGraph *si_tilt = new TGraph( vtilt.size(), &(vtilt[0]), &(vpeak[0]) ); // sim
  si_tilt->SetLineColor(2);
  si_tilt->SetLineWidth(3);
  si_tilt->SetMarkerColor(2);
  si_tilt->Draw("PL"); // without axis option: overlay

  leg->Draw();

}
