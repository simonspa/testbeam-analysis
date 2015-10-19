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
#include "plotter.C"

using namespace std;

bool cmslogo = false;

void controlplots() {
  std::cout << "Run controlplots(histogram dir)" << std::endl;
}

void controlplots(const char* inputdir, int chip, int startrun, int stoprun) {

  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.2);

  TCanvas *c1 = new TCanvas("c1","ncols",700,700);
  TProfile *ncols = new TProfile("ncols"," ",130,0,85,0,60,"");

  TCanvas *c2 = new TCanvas("c2","ncolstan",700,700);
  TProfile *ncolstan = new TProfile("ncolstan"," ",530,0,10,0,60,"");

  TCanvas *c3 = new TCanvas("c3","clustercharge",700,700);
  TProfile *clustercharge = new TProfile("clustercharge"," ",530,1,10,0,300,"");

  TCanvas *c4 = new TCanvas("c4","clustercharge",700,700);
  TProfile *clustercharge_tilt = new TProfile("clustercharge_tilt"," ",130,0,85,0,300,"");

  TCanvas *c5 = new TCanvas("c5","cmsq0f vs tilt",700,700);
  TProfile *clusterchargenorm_tilt = new TProfile("clusterchargenorm_tilt"," ",130,0,85,0,300,"");

  gStyle->SetOptStat(0);

  int nruns = 0, nclusters = 0;

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
    if(chip == 506) gDirectory->GetObject("cmsncol",h);
    else gDirectory->GetObject("cmsnrow",h);
    if(!h) continue;

    Double_t tilt = gettilt(inputdir,*run,chip);

    Double_t ncol = h->GetMean();
    Double_t landau = fitfulllang("cmsqf");
    Double_t landau_norm = fitlang("cmsq0f");

    Double_t cos = 1/TMath::Cos(TMath::Pi()*tilt/180);
    Double_t tan = TMath::Tan(TMath::Pi()*tilt/180);

    // Collect statistics:
    nclusters += h->GetEntries();

    cout << "run" << *run << " (chip" << chip << ") ncol " << ncol << " tilt " << tilt << " lanpk " << landau << " cos " << cos << " tan " << tan << " clust " << h->GetEntries() << endl;
    ncols->Fill(tilt,ncol,1);
    ncolstan->Fill(tan,ncol,1);
    clustercharge->Fill(cos,landau,1);
    clustercharge_tilt->Fill(tilt,landau,1);
    clusterchargenorm_tilt->Fill(tilt,landau_norm,1);

    nruns++;
    delete source;
  }

  cout << nruns << " runs analyzed with " << nclusters << " clusters in total." << endl;

  int thickness = 294;
  if(chip == 506) thickness = 308;

  int threshold = 170;
  if(chip == 506) threshold = 200;

  std::vector<double> vtilt = getsimulation("tilt", chip,thickness, threshold);
  std::vector<double> vtilttan = getsimulation("tilttan", chip,thickness, threshold);
  std::vector<double> vncol = getsimulation("ncol", chip,thickness, threshold);
  std::vector<double> vpath = getsimulation("path", chip,thickness, threshold);
  std::vector<double> vpeak = getsimulation("peak", chip,thickness, threshold);
  std::vector<double> vpeaknorm = getsimulation("peaknorm", chip,thickness, threshold);

  TString fileName;
  fileName.Form("chip%i-controlplots.root",chip);
  TFile * out = TFile::Open(fileName,"RECREATE");
  gDirectory->pwd();

  c1->cd();
  TLegend * leg = new TLegend();
  setLegendStyle(leg);

  setStyleAndFillLegend(ncols,"data",leg);
  if(chip == 506) ncols->SetTitle(";#alpha [#circ];columns/cluster");
  else ncols->SetTitle(";#alpha [#circ];rows/cluster");
  ncols->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
  if(chip == 506) ncols->GetYaxis()->SetRangeUser(1, 25);
  else ncols->GetYaxis()->SetRangeUser(1, 4);
  ncols->GetYaxis()->SetTitleOffset(1.2);

  ncols->Draw();
  DrawCMSLabels(nclusters,5.6,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);
  TGraph *si = new TGraph( vtilt.size(), &(vtilt[0]), &(vncol[0]) ); // sim
  si->SetLineColor(2);
  si->SetLineWidth(3);
  si->SetMarkerColor(2);
  setStyleAndFillLegend(si,"sim",leg);
  si->Draw("PL"); // without axis option: overlay
  leg->Draw();
  c1->Write();
  
  c2->cd();
  TLegend * leg2 = new TLegend();
  setLegendStyle(leg2);

  setStyleAndFillLegend(ncolstan,"data",leg2);
  if(chip == 506) ncolstan->SetTitle(";tan(#alpha);columns/cluster");
  else ncolstan->SetTitle(";tan(#alpha);rows/cluster");
  ncolstan->GetYaxis()->SetTitleOffset(1.2);
  ncolstan->GetXaxis()->SetRangeUser(vtilttan.front(), vtilttan.back());
  if(chip == 506) ncolstan->GetYaxis()->SetRangeUser(1, 25);
  else ncolstan->GetYaxis()->SetRangeUser(1, 4);
  ncolstan->Draw();
  DrawCMSLabels(nclusters,5.6,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);
  TGraph *si_tan = new TGraph( vtilttan.size(), &(vtilttan[0]), &(vncol[0]) ); // sim
  si_tan->SetLineColor(2);
  si_tan->SetLineWidth(3);
  si_tan->SetMarkerColor(2);
  si_tan->Draw("PL"); // without axis option: overlay
  setStyleAndFillLegend(si_tan,"sim",leg2);
  // From geometry:
  TF1 *fa1, *fa2;
  if(chip == 506) {
    fa1 = new TF1("fa1","1+285/150*x",vtilttan.front(),vtilttan.back()); 
    fa2 = new TF1("fa2","1+308/150*x",vtilttan.front(),vtilttan.back()); 
    leg2->AddEntry(fa1, "Geometry: 285#mum",  "l");
    leg2->AddEntry(fa2, "Geometry: 308#mum",  "l");
  }
  else {
    fa1 = new TF1("fa1","1+285/100*x",vtilttan.front(),vtilttan.back()); 
    fa2 = new TF1("fa2","1+294/100*x",vtilttan.front(),vtilttan.back()); 
    leg2->AddEntry(fa1, "Geometry: 285#mum",  "l");
    leg2->AddEntry(fa2, "Geometry: 294#mum",  "l");
  }
  fa1->SetLineWidth(3);
  fa1->SetLineStyle(7);
  fa1->SetLineColor(kBlack);
  fa2->SetLineWidth(3);
  fa2->SetLineStyle(7);
  fa2->SetLineColor(kRed+1);
  fa1->Draw("same");
  fa2->Draw("same");
  leg2->Draw();
  c2->Write();

  c3->cd();
  TLegend * leg3 = new TLegend();
  setLegendStyle(leg3);

  setStyleAndFillLegend(clustercharge,"data",leg3);
  clustercharge->SetTitle(";1/cos(#alpha);cluster charge MPV #left[ke#right]");
  clustercharge->GetYaxis()->SetTitleOffset(1.3);
  clustercharge->GetXaxis()->SetRangeUser(vpath.front(), vpath.back());
  if(chip == 506) clustercharge->GetYaxis()->SetRangeUser(20, 250);
  else clustercharge->GetYaxis()->SetRangeUser(20, 35);
  clustercharge->Draw();
  DrawCMSLabels(nclusters,5.6,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);
  TGraph *siclust = new TGraph( vpath.size(), &(vpath[0]), &(vpeak[0]) ); // sim
  siclust->SetLineColor(2);
  siclust->SetLineWidth(3);
  siclust->SetMarkerColor(2);
  siclust->Draw("PL"); // without axis option: overlay
  setStyleAndFillLegend(siclust,"sim",leg3);
  leg3->Draw();
  c3->Write();

  c4->cd();
  TLegend * leg4 = new TLegend();
  setLegendStyle(leg4);

  setStyleAndFillLegend(clustercharge_tilt,"data",leg4);
  clustercharge_tilt->SetTitle(";#alpha [#circ];cluster charge MPV #left[ke#right]");
  clustercharge_tilt->GetYaxis()->SetTitleOffset(1.3);
  clustercharge_tilt->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
  if(chip == 506) clustercharge_tilt->GetYaxis()->SetRangeUser(20, 250);
  else clustercharge_tilt->GetYaxis()->SetRangeUser(20, 35);
  clustercharge_tilt->Draw();
  DrawCMSLabels(nclusters,5.6,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);
  TGraph *si_tilt = new TGraph( vtilt.size(), &(vtilt[0]), &(vpeak[0]) ); // sim
  si_tilt->SetLineColor(2);
  si_tilt->SetLineWidth(3);
  si_tilt->SetMarkerColor(2);
  si_tilt->Draw("PL"); // without axis option: overlay
  setStyleAndFillLegend(si_tilt,"sim",leg4);
  leg4->Draw();
  c4->Write();

  c5->cd();
  TLegend * leg5 = new TLegend();
  setLegendStyle(leg5);

  setStyleAndFillLegend(clusterchargenorm_tilt,"data",leg5);
  clusterchargenorm_tilt->SetTitle(";#alpha [#circ];norm. cluster charge MPV #left[ke#right]");
  clusterchargenorm_tilt->GetYaxis()->SetTitleOffset(1.3);
  clusterchargenorm_tilt->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
  if(chip == 506) clusterchargenorm_tilt->GetYaxis()->SetRangeUser(20, 35);
  else clusterchargenorm_tilt->GetYaxis()->SetRangeUser(20, 35);
  clusterchargenorm_tilt->Draw();
  DrawCMSLabels(nclusters,5.6,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);
  TLine *line = new TLine(vtilt.front(),22,vtilt.back(),22);
  line->SetLineColor(kRed);
  line->Draw();
  TGraph *si_tiltn = new TGraph( vtilt.size(), &(vtilt[0]), &(vpeaknorm[0]) ); // sim
  si_tiltn->SetLineColor(2);
  si_tiltn->SetLineWidth(3);
  si_tiltn->SetMarkerColor(2);
  si_tiltn->Draw("PL"); // without axis option: overlay
  setStyleAndFillLegend(si_tiltn,"sim",leg5);
  leg5->Draw();
  c5->Write();
}
