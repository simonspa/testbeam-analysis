#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
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

void resolution() {
  std::cout << "Run resolution(histogram dir)" << std::endl;
}

void resolution(const char* inputdir) {
  resolution(inputdir,0,99999);
}

void resolution(const char* inputdir, int startrun, int stoprun) {

  TCanvas *c1 = new TCanvas("c1","resolution",600,600);
  TProfile *resolution = new TProfile("resolution"," ",130,5,85,0,60,"");

  TCanvas *c2 = new TCanvas("c2","resolution",600,600);
  TProfile *resolution_tel_subtracted = new TProfile("resolution_tel_subtracted"," ",130,5,85,0,60,"");

  gStyle->SetOptStat(0);

  int nruns, nfiducial, nevents;

  // Get all runs for chip 506:
  int chip = 506;

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

    // Y Resolution in fiducial volume & Landau peak:
    //Double_t res = fitep0sigma("cmsdyfctq3",-50,50);
    Double_t res = fitep0sigma("cmsdyfctq3");
    Double_t tilt = gettilt(inputdir,*run,chip);

    // Collect statistics:
    nfiducial += ((TH1D*)gDirectory->Get("cmsdyfctq3"))->GetEntries();
    nevents += ((TH1D*)gDirectory->Get("trixylk"))->GetEntries();

    // Subtract telescope track resolution:
    Double_t res_tel_subtracted = TMath::Sqrt(res*res - restel*restel);

    cout << "run" << *run << " (chip" << chip << ") res " << res << " ressub " << res_tel_subtracted << " tilt " << tilt << endl;
    resolution->Fill(tilt,res,1);
    resolution_tel_subtracted->Fill(tilt,res_tel_subtracted,1);

    nruns++;
    delete source;
  }

  cout << nruns << " runs analyzed with " << nevents << " linked clusters and " << nfiducial << " clusters in the fiducial volume in total." << endl;

  c1->cd();
  resolution->SetTitle("CMS Pixel Resolution (y) wrt to tilt angle (from: cmsdyfctq3);tilt angle [#deg];resolution y [#mum]");
  resolution->SetMarkerStyle(20);
  resolution->SetMarkerColor(2);
  resolution->Draw();

  c2->cd();
  resolution_tel_subtracted->SetTitle("CMS Pixel Resolution (y) - Tel Res. wrt to tilt angle;tilt angle [#deg];resolution y [#mum]");
  resolution_tel_subtracted->SetMarkerStyle(20);
  resolution_tel_subtracted->SetMarkerColor(2);
  resolution_tel_subtracted->Draw();

}
