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
  gStyle->SetTitleYOffset(1);

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
  int threshold = 170;
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
  if(chip == 506) ncols->SetTitle(";#alpha [#circ];columns per cluster");
  else ncols->SetTitle(";#alpha [#circ];rows per cluster");
  ncols->SetMarkerStyle(20);
  ncols->SetMarkerColor(1);
  ncols->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
  if(chip == 506) ncols->GetYaxis()->SetRangeUser(1, 25);
  else ncols->GetYaxis()->SetRangeUser(1, 4);
  ncols->Draw();
  setStyleAndFillLegend(ncols,"data",leg);
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
  if(chip == 506) ncolstan->SetTitle(";tan(#alpha);columns per cluster");
  else ncolstan->SetTitle(";tan(#alpha);rows per cluster");
  ncolstan->SetMarkerStyle(20);
  ncolstan->SetMarkerColor(1);
  ncolstan->GetXaxis()->SetRangeUser(vtilttan.front(), vtilttan.back());
  if(chip == 506) ncolstan->GetYaxis()->SetRangeUser(1, 25);
  else ncolstan->GetYaxis()->SetRangeUser(1, 4);
  ncolstan->Draw();
  setStyleAndFillLegend(ncolstan,"data",leg2);
  DrawCMSLabels(nclusters,5.6,0.045);
  if(cmslogo) DrawPrelimLabel(1,0.045);
  TGraph *si_tan = new TGraph( vtilttan.size(), &(vtilttan[0]), &(vncol[0]) ); // sim
  si_tan->SetLineColor(2);
  si_tan->SetLineWidth(3);
  si_tan->SetMarkerColor(2);
  si_tan->Draw("PL"); // without axis option: overlay
  setStyleAndFillLegend(si_tan,"sim",leg2);
  leg2->Draw();
  c2->Write();

  c3->cd();
  TLegend * leg3 = new TLegend();
  clustercharge->SetTitle(";1/cos(#alpha);peak cluster charge");
  clustercharge->SetMarkerStyle(20);
  clustercharge->SetMarkerColor(1);
  clustercharge->GetXaxis()->SetRangeUser(vpath.front(), vpath.back());
  if(chip == 506) clustercharge->GetYaxis()->SetRangeUser(20, 160);
  else clustercharge->GetYaxis()->SetRangeUser(20, 35);
  clustercharge->Draw();
  setStyleAndFillLegend(clustercharge,"data",leg3);
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
  clustercharge_tilt->SetTitle(";#alpha [#circ];peak cluster charge");
  clustercharge_tilt->SetMarkerStyle(20);
  clustercharge_tilt->SetMarkerColor(1);
  clustercharge_tilt->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
  if(chip == 506) clustercharge_tilt->GetYaxis()->SetRangeUser(20, 160);
  else clustercharge_tilt->GetYaxis()->SetRangeUser(20, 35);
  clustercharge_tilt->Draw();
  setStyleAndFillLegend(clustercharge_tilt,"data",leg4);
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
  clusterchargenorm_tilt->SetTitle(";#alpha [#circ];peak cluster chargenorm");
  clusterchargenorm_tilt->SetMarkerStyle(20);
  clusterchargenorm_tilt->SetMarkerColor(1);
  clusterchargenorm_tilt->GetXaxis()->SetRangeUser(vtilt.front(), vtilt.back());
  if(chip == 506) clusterchargenorm_tilt->GetYaxis()->SetRangeUser(20, 35);
  else clusterchargenorm_tilt->GetYaxis()->SetRangeUser(20, 35);
  clusterchargenorm_tilt->Draw();
  setStyleAndFillLegend(clusterchargenorm_tilt,"data",leg5);
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
