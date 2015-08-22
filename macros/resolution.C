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

#include "constants.h"
#include "tools.C"
#include "plotter.C"

using namespace std;

bool draw_skwcorr = false;

void resolution() {
  std::cout << "Run resolution(histogram dir)" << std::endl;
}


void resolution(const char* inputdir, int chip, int startrun, int stoprun) {

  // Set the histogram styles:
  setHHStyle(*gStyle);

  TCanvas *c1 = new TCanvas("c1","resolution",600,600);
  TProfile *resolution = new TProfile("resolution"," ",170,0,85,0,60,"");
  TProfile *resolution_corr = new TProfile("resolution_corr"," ",170,0,85,0,60,"");

  TCanvas *c2 = new TCanvas("c2","resolution",600,600);
  TProfile *resolution_tel_subtracted = new TProfile("resolution_tel_subtracted"," ",170,0,85,0,60,"");
  TProfile *resolution_corr_tel_subtracted = new TProfile("resolution_corr_tel_subtracted"," ",170,0,85,0,60,"");

  TCanvas *c3 = new TCanvas("c3","resolution",600,600);
  TProfile *resolution_vs_eta = new TProfile("resolution_vs_eta"," ",120,0.01,3,0,60,"");
  TProfile *resolution_corr_vs_eta = new TProfile("resolution_corr_vs_eta"," ",120,0.01,3,0,60,"");

  gStyle->SetOptStat(0);

  int nruns = 0, nfiducial = 0, nevents = 0;

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

    // Y Resolution in fiducial volume & Landau peak:
    //Double_t res = fitep0sigma("cmsdyfctq4",-50,50);
    Double_t res = fitep0sigma("cmsdy0fctq4d");
    //Double_t res = fittp0sigma("cmsdyfctq4d");

    //Double_t res = getRMS96("cmsdy0fctq3");
    // Skew corrected:
    Double_t res_corr = fitep0sigma("cmsdyfctq4d");
    //Double_t res_corr = getRMS96("cmsdy0fctq4",92);
    Double_t tilt = gettilt(inputdir,*run,chip);

    // Eta: flip angle from tilt to theta from beam axis:
    Double_t eta = -TMath::Log(TMath::Tan(TMath::TwoPi()*(90-tilt)/(2*360)));

    // Collect statistics:
    nfiducial += ((TH1D*)gDirectory->Get("cmsdyfctq3"))->GetEntries();
    nevents += ((TH1D*)gDirectory->Get("trixylk"))->GetEntries();

    // Subtract telescope track resolution:
    Double_t res_tel_subtracted = TMath::Sqrt(res*res - restel*restel);
    Double_t res_corr_tel_subtracted = TMath::Sqrt(res_corr*res_corr - restel*restel);

    cout << "run" << *run << " (chip" << chip << ") res " << res << " ressub " << res_tel_subtracted << " res_corr " << res_corr_tel_subtracted << " tilt " << tilt << " eta " << eta << endl;
    //cout << *run << " " << tilt << " " << res_tel_subtracted << endl;

    resolution->Fill(tilt,res,1);
    resolution_corr->Fill(tilt,res_corr,1);

    resolution_tel_subtracted->Fill(tilt,res_tel_subtracted,1);
    resolution_corr_tel_subtracted->Fill(tilt,res_corr_tel_subtracted,1);

    resolution_vs_eta->Fill(eta,res_tel_subtracted,1);
    resolution_corr_vs_eta->Fill(eta,res_corr_tel_subtracted,1);

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
  resolution->SetTitle(";tilt angle [#degree];resolution y #left[#mum#right]");
  resolution->SetMarkerStyle(20);
  resolution->SetMarkerColor(1);
  resolution->Draw();
  resolution_corr->SetMarkerStyle(20);
  resolution_corr->SetMarkerColor(15);
  if(draw_skwcorr) resolution_corr->Draw("same");
  setStyleAndFillLegend(resolution,"data",leg);
  DrawCMSLabels(nfiducial,5.6,0.045);

  c2->cd();
  resolution_tel_subtracted->SetTitle(";tilt angle [#degree];resolution y #left[#mum#right]");
  resolution_tel_subtracted->SetMarkerStyle(20);
  resolution_tel_subtracted->SetMarkerColor(1);
  resolution_tel_subtracted->Draw();
  resolution_corr_tel_subtracted->SetMarkerStyle(20);
  resolution_corr_tel_subtracted->SetMarkerColor(15);
  if(draw_skwcorr) resolution_corr_tel_subtracted->Draw("same");
  setStyleAndFillLegend(resolution_tel_subtracted,"data",leg2);
  DrawCMSLabels(nfiducial,5.6,0.045);

  c3->cd();
  resolution_vs_eta->SetTitle(";pseudo rapidity #eta;resolution y #left[#mum#right]");
  resolution_vs_eta->SetMarkerStyle(20);
  resolution_vs_eta->SetMarkerColor(1);
  resolution_vs_eta->Draw();
  resolution_corr_vs_eta->SetMarkerStyle(20);
  resolution_corr_vs_eta->SetMarkerColor(15);
  if(draw_skwcorr) resolution_corr_vs_eta->Draw("same");
  setStyleAndFillLegend(resolution_vs_eta,"data",leg3);
  DrawCMSLabels(nfiducial,5.6,0.045);
  
  int thickness = 294;
  int threshold = 170; // maybe 160 or so?
  //if(chip == 506) thickness = 308;
  std::vector<double> vtilt = getsimulation("tilt", chip,thickness,threshold);
  std::vector<double> veta = getsimulation("eta", chip,thickness,threshold);
  std::vector<double> vres = getsimulation("res", chip,thickness,threshold);
  std::vector<double> vresskw = getsimulation("resskw", chip,thickness,threshold);

  if(!vtilt.empty()) {
    c2->cd();
    TGraph *si = new TGraph( vtilt.size(), &(vtilt[0]), &(vres[0]) ); // sim
    si->SetLineColor(2);
    si->SetLineWidth(3);
    si->SetMarkerColor(2);
    resolution_tel_subtracted->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
    si->Draw("PL"); // without axis option: overlay
    setStyleAndFillLegend(si,"sim",leg2);

    TGraph *siskw = new TGraph( vtilt.size(), &(vtilt[0]), &(vresskw[0]) ); // sim
    siskw->SetLineColor(kRed-7);
    siskw->SetLineWidth(3);
    siskw->SetMarkerSize(0);
    if(draw_skwcorr) siskw->Draw("PL"); // without axis option: overlay

    c3->cd();
    TGraph *si_eta = new TGraph( veta.size(), &(veta[0]), &(vres[0]) ); // sim
    si_eta->SetLineColor(2);
    si_eta->SetLineWidth(3);
    si_eta->SetMarkerColor(2);
    resolution_vs_eta->GetXaxis()->SetRangeUser(veta.front(), veta.back());
    si_eta->Draw("PL"); // without axis option: overlay
    setStyleAndFillLegend(si_eta,"sim",leg3);

    TGraph *si_etaskw = new TGraph( veta.size(), &(veta[0]), &(vresskw[0]) ); // sim
    si_etaskw->SetLineColor(kRed-7);
    si_etaskw->SetLineWidth(3);
    si_etaskw->SetMarkerSize(0);
    if(draw_skwcorr) si_etaskw->Draw("PL"); // without axis option: overlay

  }

  c2->cd();
  leg2->Draw();
  c3->cd();
  leg3->Draw();
}
