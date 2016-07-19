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

bool cmslogo = false;

void threshold() {
  std::cout << "Run threshold(histogram dir)" << std::endl;
}


void threshold(const char* inputdir, int chip, int startrun, int stoprun) {

  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.2);

  TCanvas *c1 = new TCanvas("c1","resolution",700,700);
  TProfile *resolution = new TProfile("resolution"," ",56,1.5,4.2,3,770,"");

  TCanvas *c2 = new TCanvas("c2","resolution",700,700);
  TProfile *resolution_tel_subtracted = new TProfile("resolution_tel_subtracted"," ",56,1.5,4.2,3,7,"");

  TCanvas *c3 = new TCanvas("c3","nrows",700,700);
  TProfile *nrows = new TProfile("nrows"," ",130,1.5,4.2,0,60,"");

  TCanvas *c4 = new TCanvas("c4","lanpk",700,700);
  TProfile *lanpk = new TProfile("lanpk"," ",130,1.5,4.2,0,60,"");

  gStyle->SetOptStat(0);

  int nruns = 0, nfiducial = 0, nevents = 0;
  Double_t tilt;
  Double_t tmp_ke, tmp_res, tmp_ressub, tmp_ncol, tmp_lanpk;

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
    if(chip == 506) gDirectory->GetObject("cmsncol",h);
    else gDirectory->GetObject("cmsnrow",h);
    if(!h) continue;

    // Y Resolution in fiducial volume & Landau peak:
    //Double_t res = fitep0sigma("cmsdyfctq4",-50,50);
    Double_t res = fitep0sigma("cmsdy0fctq4d");
    //Double_t res = fittp0sigma("cmsdyfctq4d");

    tilt = gettilt(inputdir,*run,chip,"-threshold");
    Double_t ncol = h->GetMean();

    // Collect statistics:
    nfiducial += ((TH1D*)gDirectory->Get("cmsdyfctq3"))->GetEntries();
    nevents += ((TH1D*)gDirectory->Get("trixylk"))->GetEntries();

    int threshold = gettrim(inputdir,*run,chip,"-threshold");

    // From VCal to electrons: 50e/VCal DAC
    Double_t ke = threshold*50.0/1000;

    Double_t landau = fitfulllang("cmsqf");

    // Subtract telescope track resolution:
    Double_t dz = getdz(inputdir,*run,chip,"-threshold");
    Double_t tel_resolution = getTelRes(dz);

    Double_t res_tel_subtracted = TMath::Sqrt(res*res - tel_resolution*tel_resolution);

    cout << "run" << *run << " (chip" << chip << ") res " << res << " ressub " << res_tel_subtracted << " tilt " << tilt << " electrons " << ke << " ncol " << ncol << " threshold " << threshold << endl;
    //cout << *run << " " << tilt << " " << res_tel_subtracted << endl;
    //cout << " -> dz = " << dz << ", sigma_tel = " << tel_resolution << endl;

    nrows->Fill(ke,ncol);
    lanpk->Fill(ke,landau);
    resolution->Fill(ke,res,1);
    resolution_tel_subtracted->Fill(ke,res_tel_subtracted,1);
    tmp_ke = ke; tmp_res = res; tmp_ressub = res_tel_subtracted; tmp_ncol = ncol; tmp_lanpk = landau;
    nruns++;
    delete source;
  }

  cout << nruns << " runs analyzed with " << nevents << " linked clusters and " << nfiducial << " clusters in the fiducial volume in total." << endl;

  resolution->Fill(tmp_ke,tmp_res+0.1,1);
  resolution_tel_subtracted->Fill(tmp_ke,tmp_ressub+0.2,1);
  nrows->Fill(tmp_ke,tmp_ncol+0.01,1);
  lanpk->Fill(tmp_ke,tmp_lanpk+0.01,1);


  tilt++;
  int thickness = 300;
  // Read simulations:
  std::vector<double> vthr = getthresholds("thr", tilt, thickness, chip);
  std::vector<double> vres = getthresholds("res", tilt, thickness, chip);
  std::vector<double> vreserr_up = getthresholds("res", tilt, 308, chip);
  std::vector<double> vreserr_down = getthresholds("res", tilt, 294, chip);
  std::vector<double> vncol = getthresholds("ncol", tilt, thickness, chip);
  std::vector<double> vncolerr_up = getthresholds("ncol", tilt, 308, chip);
  std::vector<double> vncolerr_down = getthresholds("ncol", tilt, 294, chip);
  std::vector<double> vlanpk = getthresholds("lanpk", tilt, thickness, chip);

  std::vector<double> nullvec;
  std::vector<double> vreserr;
  std::vector<double> vncolerr;

  cout << "vres " << vres.size() << " vreserr down " << vreserr_down.size() << " vreserr up " << vreserr_up.size() << endl;
  int min = std::min(vres.size(),std::min(vreserr_up.size(),vreserr_down.size()));

  for(int i = 0; i < min; i++) {
    double err_thick = max(fabs(vreserr_down.at(i)-vres.at(i)),fabs(vreserr_up.at(i)-vres.at(i)));
    vreserr.push_back(err_thick);
    cout << vreserr.at(i) << " = " 
	 << max(fabs(vreserr_down.at(i)-vres.at(i)),fabs(vreserr_up.at(i)-vres.at(i))) 
	 << endl;
    nullvec.push_back(0);
  }

  min = std::min(vncol.size(),std::min(vncolerr_up.size(),vncolerr_down.size()));
  for(int i = 0; i < min; i++) {
    double err_thick = max(fabs(vncolerr_down.at(i)-vncol.at(i)),fabs(vncolerr_up.at(i)-vncol.at(i)));
    vncolerr.push_back(err_thick);
    cout << vncolerr.at(i) << " = " 
	 << max(fabs(vncolerr_down.at(i)-vncol.at(i)),fabs(vncolerr_up.at(i)-vncol.at(i))) 
	 << endl;
    nullvec.push_back(0);
  }


  TString fileName;
  fileName.Form("chip%i-threshold.root",chip);
  TFile * out = TFile::Open(fileName,"RECREATE");
  gDirectory->pwd();

  TLegend *leg = new TLegend();
  TLegend *leg2 = new TLegend();
  TLegend *leg3 = new TLegend();
  TLegend *leg4 = new TLegend();
  setLegendStyle(leg);
  setLegendStyle(leg2);
  setLegendStyle(leg3);
  setLegendStyle(leg4);

  c1->cd();
  resolution->SetTitle(";pixel threshold [ke];resolution x #left[#mum#right]");
  resolution->SetMarkerStyle(20);
  resolution->SetMarkerColor(1);
  resolution->Draw("e");
  setStyleAndFillLegend(resolution,"data",leg);
  DrawCMSLabels(nfiducial,5.2,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);
  c1->Write();

  c2->cd();
  resolution_tel_subtracted->SetTitle(";pixel threshold [ke];resolution x #left[#mum#right]");
  resolution_tel_subtracted->SetMarkerStyle(20);
  resolution_tel_subtracted->SetMarkerColor(1);
  resolution_tel_subtracted->GetYaxis()->SetRangeUser(0, 15);
  resolution_tel_subtracted->Draw("e");
  setStyleAndFillLegend(resolution_tel_subtracted,"data",leg2);
  DrawCMSLabels(nfiducial,5.2,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);

  c3->cd();
  if(chip == 506) nrows->SetTitle(";pixel threshold [ke];<pixels/cluster>");
  else nrows->SetTitle(";pixel threshold [ke];#LTpixels/cluster#GT");
  nrows->SetMarkerStyle(20);
  nrows->SetMarkerColor(1);
  nrows->GetXaxis()->SetRangeUser(vthr.front(), vthr.back());
  if(chip == 506) nrows->GetYaxis()->SetRangeUser(1, 21);
  else nrows->GetYaxis()->SetRangeUser(1, 4);
  nrows->Draw("e");
  setStyleAndFillLegend(nrows,"data",leg3);
  DrawCMSLabels(nfiducial,5.2,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);

  c4->cd();
  lanpk->SetTitle(";pixel threshold [ke];Landau MPV #left[ke#right]");
  lanpk->SetMarkerStyle(20);
  lanpk->SetMarkerColor(1);
  lanpk->Draw("e");
  lanpk->GetYaxis()->SetRangeUser(20, 26);
  setStyleAndFillLegend(lanpk,"data",leg4);
  DrawCMSLabels(nfiducial,5.2,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);

  if(!vthr.empty()) {
    c2->cd();
    TGraphAsymmErrors *si = new TGraphAsymmErrors( vthr.size(), &(vthr[0]), &(vres[0]), &(nullvec[0]), &(nullvec[0]), &(vreserr[0]), &(vreserr[0]) ); // sim
    si->SetLineColor(2);
    si->SetLineWidth(3);
    si->SetMarkerColor(2);
    resolution_tel_subtracted->GetXaxis()->SetRangeUser(vthr.front(), vthr.back());
    setStyleAndFillLegend(si,"sim",leg2);
    si->Draw("PC4"); // without axis option: overlay
  }
  c2->cd();
  leg2->Draw();
  c2->Write();

  if(!vthr.empty()) {
    c3->cd();
    TGraphAsymmErrors *si = new TGraphAsymmErrors( vthr.size(), &(vthr[0]), &(vncol[0]), &(nullvec[0]), &(nullvec[0]), &(vncolerr[0]), &(vncolerr[0]) ); // sim
    si->SetLineColor(2);
    si->SetLineWidth(3);
    si->SetMarkerColor(2);
    setStyleAndFillLegend(si,"sim",leg3);
    si->Draw("PC4"); // without axis option: overlay
  }
  c3->cd();
  leg3->Draw();
  c3->Write();

  if(!vthr.empty()) {
    c4->cd();
    TGraph *si = new TGraph( vthr.size(), &(vthr[0]), &(vlanpk[0]) ); // sim
    si->SetLineColor(2);
    si->SetLineWidth(3);
    si->SetMarkerColor(2);
    setStyleAndFillLegend(si,"sim",leg4);
    si->Draw("PL"); // without axis option: overlay
  }
  c4->cd();
  leg4->Draw();
  c4->Write();

}
