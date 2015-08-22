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

void efficiency() {
  std::cout << "Run efficiency(histogram dir)" << std::endl;
}

void efficiency(const char* inputFilesDirectory) {

  TCanvas *effmap = new TCanvas("c1","efficiency_map",600,600);
  TProfile2D * effvsxy = new TProfile2D("effvsxy","eff", 60, -4.5, 4.5, 90, -4.5, 4.5,);

  gStyle->SetOptStat(0);

  // Get all runs for chip 500:
  std::vector<int> runs = getruns(500);
  for(std::vector<int>::iterator run = runs.begin(); run != runs.end(); run++) {

      TString fileName;
      fileName += inputFilesDirectory;
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
      Double_t res = fitgp0("cmsdyfctq3d");
      nfiducial += ((TH1D*)gDirectory->Get("cmsdyfctq3d"))->GetEntries();
      // Number of Rows per Cluster:
      Double_t nrow = ((TH1D*)gDirectory->Get("cmsnrowq"))->GetMean();
      // Landau Peak Position
      Double_t land = fitlang("cmsqf");

      nevents += ((TH1D*)gDirectory->Get("trixylk"))->GetEntries();

      cout << "run" << *run << " (chip" << chip << ") res " << res << " tilt " << tilt << " land " << land << " nrow " << nrow << endl;

      if(chip == 500) hres500->Fill(tilt,res,1);
      else hres504->Fill(tilt,res,1);

      if(chip == 500) hnrow500->Fill(tilt,nrow,1);
      else hnrow504->Fill(tilt,nrow,1);

      if(chip == 500) hland500->Fill(tilt,land,1);
      else hland504->Fill(tilt,land,1);

      Double_t cos_tilt = (1/cos(tilt*2*TMath::Pi()/360));
      if(chip == 500) hland2_500->Fill(cos_tilt,land,1);
      else hland2_504->Fill(cos_tilt,land,1);

      nruns++;
      delete source;
    }
  }

  cout << nruns << " runs analyzed with " << nevents << " linked clusters and " << nfiducial << " clusters in the fiducial volume in total." << endl;

  c1->cd();
  hres500->SetTitle("CMS Pixel Resolution (y) wrt to tilt angle (from: cmsdyfctq3d);tilt angle [#deg];resolution y [#mum]");
  hres500->SetMarkerStyle(20);
  hres500->SetMarkerColor(2);
  hres500->Draw();
  hres504->SetMarkerStyle(22);
  hres504->SetMarkerColor(4);
  hres504->Draw("same");

  c2->cd();
  hnrow500->SetTitle("CMS Pixel No. of Rows wrt to tilt angle (from: cmsnrowq);tilt angle [#deg];row per cluster");
  hnrow500->SetMarkerStyle(20);
  hnrow500->SetMarkerColor(2);
  hnrow500->Draw();
  hnrow504->SetMarkerStyle(22);
  hnrow504->SetMarkerColor(4);
  hnrow504->Draw("same");

  c3->cd();
  hland504->SetTitle("CMS Pixel Landau Pos. (from: cmsqf);tilt angle [#deg];landau pos [ke]");
  hland504->SetMarkerStyle(22);
  hland504->SetMarkerColor(4);
  hland504->Draw();
  hland500->SetMarkerStyle(20);
  hland500->SetMarkerColor(2);
  hland500->Draw("same");


  c4->cd();
  hland2_504->SetTitle("CMS Pixel Landau Pos. (from: cmsqf);tilt: 1/cos(#phi);landau pos [ke]");
  hland2_504->SetMarkerStyle(22);
  hland2_504->SetMarkerColor(4);
  hland2_504->Draw();
  hland2_500->SetMarkerStyle(20);
  hland2_500->SetMarkerColor(2);
  hland2_500->Draw("same");


  in.close();
  MyFile->Write();
  MyFile->Close();
  //  getchar();
}
