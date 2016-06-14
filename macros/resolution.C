#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
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

bool cmslogo = false;

void resolution() {
  std::cout << "Run resolution(histogram dir)" << std::endl;
}


void resolution(const char* inputdir, int chip, int startrun, int stoprun, bool draw_skwcorr = true, string postfix = "", int threshold = 170) {

  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.1);

  TCanvas *c1 = new TCanvas("c1","resolution",700,700);
  TProfile *resolution = new TProfile("resolution"," ",170,0,85,0,60,"");
  TProfile *resolution_corr = new TProfile("resolution_corr"," ",170,0,85,0,60,"");

  TCanvas *c2 = new TCanvas("c2","resolution",700,700);
  TProfile *resolution_tel_subtracted = new TProfile("resolution_tel_subtracted"," ",170,0,85,0,60,"");
  TProfile *resolution_corr_tel_subtracted = new TProfile("resolution_corr_tel_subtracted"," ",170,0,85,0,60,"");

  TCanvas *c3 = new TCanvas("c3","resolution",700,700);
  TProfile *resolution_vs_eta = new TProfile("resolution_vs_eta"," ",120,0.01,3,0,60,"");
  TProfile *resolution_corr_vs_eta = new TProfile("resolution_corr_vs_eta"," ",120,0.01,3,0,60,"");

  TCanvas *c4 = new TCanvas("c4","resolution_improvement",700,700);
  TProfile *resolution_improvement = new TProfile("resolution_improvement"," ",170,0,85,-10,100,"");

  gStyle->SetOptStat(0);

  int nruns = 0, nfiducial = 0, nevents = 0;

  // Get all runs for given chip:
  std::vector<int> runs = getruns(inputdir,chip,postfix);
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
    Double_t res = fitep0sigma("cmsdy0fctq4d");
    //Double_t res = getRMS96("cmsdy0fctq3");

    // Skew corrected:
    Double_t res_corr = fitep0sigma("cmsdyfctq4d");
    //Double_t res_corr = getRMS96("cmsdy0fctq4",92);

    Double_t tilt = gettilt(inputdir,*run,chip,postfix);

    // Eta: flip angle from tilt to theta from beam axis:
    Double_t eta = -TMath::Log(TMath::Tan(TMath::TwoPi()*(90-tilt)/(2*360)));

    // Collect statistics:
    nfiducial += ((TH1D*)gDirectory->Get("cmsdyfctq3"))->GetEntries();
    nevents += ((TH1D*)gDirectory->Get("trixylk"))->GetEntries();

    // Subtract telescope track resolution:
    Double_t dz = getdz(inputdir,*run,chip,postfix);
    Double_t tel_resolution = getTelRes(dz);

    Double_t res_tel_subtracted = TMath::Sqrt(res*res - tel_resolution*tel_resolution);
    Double_t res_corr_tel_subtracted = TMath::Sqrt(res_corr*res_corr - tel_resolution*tel_resolution);

    Double_t impr = (1-res_corr_tel_subtracted/res_tel_subtracted)*100;

    cout << "run" << *run << " (chip" << chip << ") res " << res 
	 << " ressub " << res_tel_subtracted 
	 << " res_corr " << res_corr_tel_subtracted 
	 << " (" << impr << "%)"
	 << " tilt " << tilt 
	 << " eta " << eta << endl;

    resolution->Fill(tilt,res,1);
    resolution_corr->Fill(tilt,res_corr,1);

    resolution_tel_subtracted->Fill(tilt,res_tel_subtracted,1);
    resolution_corr_tel_subtracted->Fill(tilt,res_corr_tel_subtracted,1);
    resolution_improvement->Fill(tilt,impr,1);

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

  TString fileName;
  fileName.Form("chip%i-resolution.root",chip);
  TFile * out = TFile::Open(fileName,"RECREATE");
  gDirectory->pwd();

  c1->cd();
  if(chip == 506) resolution->SetTitle(";tilt angle #omega [#circ];resolution y #left[#mum#right]");
  else resolution->SetTitle(";tilt angle #alpha [#circ];resolution x #left[#mum#right]");
  resolution->SetMarkerStyle(20);
  resolution_corr->SetMarkerStyle(20);
  if(draw_skwcorr) {
    resolution->SetMarkerColor(kGray);
    resolution_corr->SetMarkerColor(kBlack);
  }
  else { resolution->SetMarkerColor(kBlack); }

  resolution->Draw();
  if(draw_skwcorr) {
    resolution_corr->Draw("same");
    setStyle(resolution_corr,"data");
    leg->AddEntry(resolution, "Data", "p");
    leg->AddEntry(resolution_corr, "Data (corrected)", "p");
  }
  else {
    setStyleAndFillLegend(resolution,"data",leg);
  }
  DrawCMSLabels(nfiducial,5.6,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);
  c1->Write();

  c2->cd();
  if(chip == 506) resolution_tel_subtracted->SetTitle(";tilt angle #omega [#circ];resolution y #left[#mum#right]");
  else resolution_tel_subtracted->SetTitle(";tilt angle #alpha [#circ];resolution x #left[#mum#right]");
  resolution_tel_subtracted->SetMarkerStyle(20);
  resolution_corr_tel_subtracted->SetMarkerStyle(20);
  if(draw_skwcorr) {
    resolution_tel_subtracted->SetMarkerColor(kGray);
    resolution_corr_tel_subtracted->SetMarkerColor(kBlack);
  }
  else { resolution_tel_subtracted->SetMarkerColor(kBlack); }

  resolution_tel_subtracted->Draw();
  if(draw_skwcorr) {
    resolution_corr_tel_subtracted->Draw("same");
    setStyle(resolution_corr_tel_subtracted,"data");
    leg2->AddEntry(resolution_tel_subtracted, "Data", "p");
    leg2->AddEntry(resolution_corr_tel_subtracted, "Data (corrected)", "p");
  }
  else {
    setStyleAndFillLegend(resolution_tel_subtracted,"data",leg2);
  }
  DrawCMSLabels(nfiducial,5.6,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);

  if(chip == 506) {
    c3->cd();
    resolution_vs_eta->SetTitle(";pseudo rapidity #eta;resolution y #left[#mum#right]");
    resolution_vs_eta->SetMarkerStyle(20);
    resolution_corr_vs_eta->SetMarkerStyle(20);
    if(draw_skwcorr) {
      resolution_vs_eta->SetMarkerColor(kGray);
      resolution_corr_vs_eta->SetMarkerColor(kBlack);
    }
    else { resolution_vs_eta->SetMarkerColor(kBlack); }

    resolution_vs_eta->Draw();
    if(draw_skwcorr) {
      resolution_corr_vs_eta->Draw("same");
      setStyle(resolution_corr_vs_eta,"data");
      leg3->AddEntry(resolution_vs_eta, "Data", "p");
      leg3->AddEntry(resolution_corr_vs_eta, "Data (corrected)", "p");
    }
    else {
      setStyleAndFillLegend(resolution_vs_eta,"data",leg3);
    }
    DrawCMSLabels(nfiducial,5.6,0.045);
    if(cmslogo) DrawPrelimLabel(1,0.045);
  }


  // Simulations:  
  int thickness = 300;
  //if(chip == 506) thickness = 308;

  //int threshold = 400;//170;
  if(chip == 506) threshold = 200;

  std::vector<double> vtilt = getsimulation("tilt", chip,thickness,threshold);
  std::vector<double> veta = getsimulation("eta", chip,thickness,threshold);
  std::vector<double> vres = getsimulation("res", chip,thickness,threshold);
  std::vector<double> vreserr_down = getsimulation("res", chip,294,threshold);
  std::vector<double> vreserr_up = getsimulation("res", chip,308,threshold);
  std::vector<double> vresthr_down = getsimulation("res", chip,thickness,threshold-30);
  std::vector<double> vresthr_up = getsimulation("res", chip,thickness,threshold+50);
  std::vector<double> vresskw = getsimulation("resskw", chip,thickness,threshold);

  std::vector<double> nullvec;
  std::vector<double> vreserr;

  cout << "vres " << vres.size() << " vreserr down " << vreserr_down.size() << " vreserr up " << vreserr_up.size() << endl;
  cout << "vres " << vres.size() << " vresthr down " << vresthr_down.size() << " vresthr up " << vresthr_up.size() << endl;
  int min = std::min(vres.size(),std::min(std::min(vreserr_up.size(),vreserr_down.size()),std::min(vresthr_up.size(),vresthr_down.size())));

  for(int i = 0; i < min; i++) {
    double err_thick = max(fabs(vreserr_down.at(i)-vres.at(i)),fabs(vreserr_up.at(i)-vres.at(i)));
    double err_thresh = max(fabs(vresthr_down.at(i)-vres.at(i)),fabs(vresthr_up.at(i)-vres.at(i)));
    double err_total = sqrt(err_thick*err_thick + err_thresh*err_thresh);
    vreserr.push_back(err_total);
    cout << vreserr.at(i) << " = sqrt(" 
	 << max(fabs(vreserr_down.at(i)-vres.at(i)),fabs(vreserr_up.at(i)-vres.at(i))) 
	 << "**2 + " 
	 << max(fabs(vresthr_down.at(i)-vres.at(i)),fabs(vresthr_up.at(i)-vres.at(i))) 
	 << "**2)" <<endl;
    nullvec.push_back(0);
  }
  cout << endl;

  if(!vtilt.empty()) {
    c2->cd();
    TGraphAsymmErrors *si = new TGraphAsymmErrors( vtilt.size(), &(vtilt[0]), &(vres[0]), &(nullvec[0]), &(nullvec[0]), &(vreserr[0]), &(vreserr[0]) ); // sim
    TGraph *siskw = new TGraph( vtilt.size(), &(vtilt[0]), &(vresskw[0]) ); // sim
    si->SetLineWidth(3);
    siskw->SetLineWidth(3);
    si->SetMarkerSize(0);
    siskw->SetMarkerSize(0);

    if(draw_skwcorr) {
      si->SetLineColor(kGray+1);
      siskw->SetLineColor(2);
    }
    else { si->SetLineColor(2); }

    resolution_tel_subtracted->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());

    si->Draw("PL3"); // without axis option: overlay
    resolution_tel_subtracted->Draw("same");

    if(draw_skwcorr) {
      siskw->Draw("PL"); // without axis option: overlay
      setStyle(siskw,"sim");
      leg2->AddEntry(si, "Simulation", "l");
      leg2->AddEntry(siskw, "Simulation (corrected)", "l");
    }
    else { 
      setStyleAndFillLegend(si,"sim",leg2); 
      si->SetFillStyle(3005);
      si->SetFillColor(kRed+1);
    }

    if(chip == 506) {
      c3->cd();
      TGraphAsymmErrors *si_eta = new TGraphAsymmErrors( veta.size(), &(veta[0]), &(vres[0]), &(nullvec[0]), &(nullvec[0]), &(vreserr[0]), &(vreserr[0]) ); // sim
      TGraph *si_etaskw = new TGraph( veta.size(), &(veta[0]), &(vresskw[0]) ); // sim
      si_eta->SetLineWidth(3);
      si_eta->SetMarkerSize(0);
      si_etaskw->SetLineWidth(3);
      si_etaskw->SetMarkerSize(0);

      if(draw_skwcorr) {
	si_eta->SetLineColor(kGray+1);
	si_etaskw->SetLineColor(2);
      }
      else { si_eta->SetLineColor(2); }

      resolution_vs_eta->GetXaxis()->SetRangeUser(veta.front(), veta.back());
      si_eta->Draw("PL3"); // without axis option: overlay
      resolution_vs_eta->Draw("same");

      if(draw_skwcorr) {
	si_etaskw->Draw("PL"); // without axis option: overlay
	setStyle(si_etaskw,"sim");
	leg3->AddEntry(si_eta, "Simulation", "l");
	leg3->AddEntry(si_etaskw, "Simulation (corrected)", "l");
      }
      else { 
	setStyleAndFillLegend(si_eta,"sim",leg3); 
	si_eta->SetFillStyle(3005);
	si_eta->SetFillColor(kRed+1);
      }
    }
  }

  c2->cd();
  leg2->Draw();
  c2->Write();

  if(chip == 506) {
    c3->cd();
    leg3->Draw();
    c3->Write();
  }

  if(draw_skwcorr) {
    c4->cd();
    if(chip == 506) resolution_improvement->SetTitle(";tilt angle #omega [#circ];rel. improvement #left[%#right]");
    else resolution_improvement->SetTitle(";tilt angle #alpha [#circ];rel. improvement #left[%#right]");
    setStyle(resolution_improvement,"data");

    //resolution_improvement->SetMarkerStyle(0);
    resolution_improvement->SetLineColor(kGray+2);
    //resolution_improvement->SetLineWidth(0);
    resolution_improvement->SetFillStyle(3004);
    resolution_improvement->SetFillColor(kBlack);
    if(!vtilt.empty()) resolution_improvement->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
    resolution_improvement->Draw();
    DrawCMSLabels(nfiducial,5.6,0.045);
    if(cmslogo) DrawPrelimLabel(1,0.045);
  }
  c4->Write();

}
