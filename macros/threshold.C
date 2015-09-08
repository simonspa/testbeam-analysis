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
#include "TString.h"
#include "TROOT.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include "constants.h"
#include "tools.C"
#include "plotter.C"

using namespace std;

bool draw_skwcorr = false;

void threshold() {
  std::cout << "Run threshold(histogram dir)" << std::endl;
}


void threshold(const char* inputdir, int chip, int startrun, int stoprun) {

  // Set the histogram styles:
  setHHStyle(*gStyle);

  TCanvas *c1 = new TCanvas("c1","resolution",600,600);
  TProfile *resolution = new TProfile("resolution"," ",56,1500,4000,3,770,"");
  TProfile *resolution_corr = new TProfile("resolution_corr"," ",56,1500,4000,3,7,"");

  TCanvas *c2 = new TCanvas("c2","resolution",600,600);
  TProfile *resolution_tel_subtracted = new TProfile("resolution_tel_subtracted"," ",56,1500,4000,3,7,"");
  TProfile *resolution_corr_tel_subtracted = new TProfile("resolution_corr_tel_subtracted"," ",56,1500,4000,3,7,"");

  gStyle->SetOptStat(0);

  int nruns = 0, nfiducial = 0, nevents = 0;
  Double_t tilt;

  // Get all runs for given chip:
  std::vector<int> runs = getruns(inputdir,chip,"-threshold");
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
    //Double_t res = fitep0sigma("cmsdyfctq4",-50,50);
    Double_t res = fitep0sigma("cmsdy0fctq4d");
    //Double_t res = fittp0sigma("cmsdyfctq4d");

    //Double_t res = getRMS96("cmsdy0fctq3");
    // Skew corrected:
    Double_t res_corr = fitep0sigma("cmsdyfctq4d");
    //Double_t res_corr = getRMS96("cmsdy0fctq4",92);
    tilt = gettilt(inputdir,*run,chip,"-threshold");

    // Collect statistics:
    nfiducial += ((TH1D*)gDirectory->Get("cmsdyfctq3"))->GetEntries();
    nevents += ((TH1D*)gDirectory->Get("trixylk"))->GetEntries();

    int threshold = gettrim(inputdir,*run,chip,"-threshold");

    // From VCal to electrons: 50e/VCal DAC
    Double_t ke = threshold*50;

    // Subtract telescope track resolution:
    Double_t dz = getdz(inputdir,*run,chip,"-threshold");
    Double_t tel_resolution = getTelRes(dz);

    Double_t res_tel_subtracted = TMath::Sqrt(res*res - tel_resolution*tel_resolution);
    Double_t res_corr_tel_subtracted = TMath::Sqrt(res_corr*res_corr - tel_resolution*tel_resolution);

    cout << "run" << *run << " (chip" << chip << ") res " << res << " ressub " << res_tel_subtracted << " res_corr " << res_corr_tel_subtracted << " tilt " << tilt << " electrons " << ke << " threshold " << threshold << endl;
    //cout << *run << " " << tilt << " " << res_tel_subtracted << endl;
    //cout << " -> dz = " << dz << ", sigma_tel = " << tel_resolution << endl;

    resolution->Fill(ke,res,1);
    resolution_corr->Fill(ke,res_corr,1);

    resolution_tel_subtracted->Fill(ke,res_tel_subtracted,1);
    resolution_corr_tel_subtracted->Fill(ke,res_corr_tel_subtracted,1);

    nruns++;
    delete source;
  }

  cout << nruns << " runs analyzed with " << nevents << " linked clusters and " << nfiducial << " clusters in the fiducial volume in total." << endl;

  TLegend *leg = new TLegend();
  TLegend *leg2 = new TLegend();
  TLegend *leg3 = new TLegend();
  setLegendStyle(leg);
  setLegendStyle(leg2);
  setLegendStyle(leg3);

  c1->cd();
  resolution->SetTitle(";threshold [e];resolution x #left[#mum#right]");
  resolution->SetMarkerStyle(20);
  resolution->SetMarkerColor(1);
  resolution->Draw();
  resolution_corr->SetMarkerStyle(20);
  resolution_corr->SetMarkerColor(15);
  if(draw_skwcorr) resolution_corr->Draw("same");
  setStyleAndFillLegend(resolution,"data",leg);
  DrawCMSLabels(nfiducial,5.2,0.045);

  c2->cd();
  resolution_tel_subtracted->SetTitle(";threshold [e];resolution x #left[#mum#right]");
  resolution_tel_subtracted->SetMarkerStyle(20);
  resolution_tel_subtracted->SetMarkerColor(1);
  resolution_tel_subtracted->Draw();
  resolution_corr_tel_subtracted->SetMarkerStyle(20);
  resolution_corr_tel_subtracted->SetMarkerColor(15);
  if(draw_skwcorr) resolution_corr_tel_subtracted->Draw("same");
  setStyleAndFillLegend(resolution_tel_subtracted,"data",leg2);
  DrawCMSLabels(nfiducial,5.2,0.045);


  int thickness = 294;

  std::vector<double> vthr;
  std::vector<double> vres;

  // Read thresholds
  for(Int_t thr = 150; thr < 420; thr += 10) {

    TString myfile;
    myfile.Form("simulation/sim%i_%iskw_thr%i.dat",thickness,chip,thr);
  
    ifstream SIMstream( myfile );
    if( !SIMstream ) { continue;}

    // Read file by lines:
    string rl;
    double tlt;
    int run;
    int nev;
    double ry;
    double ry_skwcorr;
    double ncol;
    double lanpk;
    double turn;
    double edge;

    while( SIMstream.good() && ! SIMstream.eof() ) {

      getline( SIMstream, rl ); // read one line  = event into string
      istringstream simrun( rl ); // tokenize string

      simrun >> run;
      simrun >> tlt; // [deg]
      simrun >> turn;
      simrun >> nev;
      simrun >> ry;
      simrun >> ry_skwcorr;
      simrun >> ncol;
      simrun >> lanpk;
      simrun >> edge;

      if(tilt < (tlt+0.5)) {
	vthr.push_back(10*thr);
	vres.push_back(sqrt( ry*ry - restel_sim*restel_sim )); // subtract telescope
	cout << "sim tilt " << tlt << " res " << ry << " ressub " << sqrt( ry*ry - restel_sim*restel_sim ) << " electrons " << (thr*10) << endl;
	break;
      }
    } // while lines
  }

  if(!vthr.empty()) {
    c2->cd();
    TGraph *si = new TGraph( vthr.size(), &(vthr[0]), &(vres[0]) ); // sim
    si->SetLineColor(2);
    si->SetLineWidth(3);
    si->SetMarkerColor(2);
    resolution_tel_subtracted->GetXaxis()->SetRangeUser(vthr.front(), vthr.back());
    si->Draw("PL"); // without axis option: overlay
    setStyleAndFillLegend(si,"sim",leg2);
  }

  c2->cd();
  leg2->Draw();

}
