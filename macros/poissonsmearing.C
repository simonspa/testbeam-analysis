#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TRandom3.h"
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

#include "constants.h"
#include "tools.C"
#include "plotter.C"

using namespace std;

void poissonsmearing() {
  std::cout << "Run poissonsmearing(histogram dir)" << std::endl;
}


void poissonsmearing(const char* inputdir, int myrun, char* histogram, int nexperiments = 1000) {

  // Get random number generator:
  TRandom3 *myrnd = new TRandom3();

  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.1);
  gStyle->SetOptStat(0);


  TCanvas *c1 = new TCanvas("c1","poissonsmearing",700,700);
  //TH1D *sigmas = new TH1D("sigmas"," ",200,0.9*orig_res,1.1*orig_res);
  TH1D *sigmas = new TH1D("sigmas","",200,6,7.5);

  TString fileName;
  fileName += inputdir;
  if( !fileName.EndsWith("/") ) fileName += "/";
  fileName += "run" + ZeroPadNumber(myrun,6) + "-analysis.root";

  TFile *source;
  if (!gSystem->AccessPathName( fileName )) source = TFile::Open(fileName);
  else return;

  // get the original histogram:
  source->cd("MyEUTelAnalysisCMSPixel");
  TH1D *original;
  gDirectory->GetObject(histogram,original);
  if(!original) return;

  // Fit original for reference:
  Double_t orig_res = fitep0sigma(histogram);
  cout << "orig: sigma = " << orig_res << ", histo " << 0.9*orig_res << "-" << 1.1*orig_res << endl;
  // New canvas in +-10% around original:

  // Repeat this pseudo experiment 1k times:
  for(size_t i = 0; i < nexperiments; i++) {

    // new histogram, clone of the original:
    TH1D * smeared = (TH1D*)original->Clone("smeared");

    // loop over the bins:
    for(Int_t bin = 1; bin <= original->GetNbinsX(); bin++) {
      // Get random number from Poisson distribution with mean of original bin content:
      Double_t content = myrnd->Poisson(original->GetBinContent(bin));
      //cout << " " << i << " bin " << bin << " orig " << original->GetBinContent(bin) << " rnd " << content << endl;
      smeared->SetBinContent(bin,content);
    }

    Double_t res = fitep0sigma("smeared");

    if(i%1000 == 0) cout << "i" << i << " sigma = " << res << endl;
    sigmas->Fill(res,1);

    delete smeared;
  }

  delete source;

  cout << nexperiments << " pseudo experiments performed." << endl;

  c1->cd();
  cout << nexperiments << " pseudo experiments performed." << endl;
  sigmas->SetTitle(";#sigma [#mum];entries");
  sigmas->GetXaxis()->SetLabelFont(42);
  sigmas->GetYaxis()->SetLabelFont(42);
  sigmas->GetXaxis()->SetTitleFont(42);
  sigmas->GetYaxis()->SetTitleFont(42);
  sigmas->GetXaxis()->SetTitleSize(0.05);
  sigmas->GetYaxis()->SetTitleSize(0.05);
  sigmas->GetXaxis()->SetTitleOffset(1.08);
  sigmas->GetYaxis()->SetTitleOffset(1.2);
  sigmas->GetXaxis()->SetLabelOffset(0.007);
  sigmas->GetYaxis()->SetLabelOffset(0.007);

  sigmas->SetMarkerStyle(20);
  sigmas->SetMarkerSize(1.2);
  sigmas->SetLineWidth(3);
  sigmas->SetLineColor(kBlack);
  sigmas->SetMarkerColor(kBlack);
  sigmas->SetFillStyle(3005);
  sigmas->SetFillColor(kBlack);
  sigmas->Draw();

  TLegend * leg =  new TLegend();
  DrawFreeCMSLabels(Form("%2.1f #times 10^{3} pseudo exp.",1.0*nexperiments/1000),0.045);

  sigmas->Fit("gaus");
  TF1 *fit = sigmas->GetFunction("gaus");
  fit->SetLineWidth(3);
  fit->SetLineColor(kRed+1);
  cout << "statistical uncertainty: " << fit->GetParameter(2) << endl;
  leg->Draw();

}
