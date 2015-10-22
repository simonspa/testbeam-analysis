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

void fitskw_data() {
  cout << "inputdir!" << endl;
}

void fitskw_data(TString inputdir, int chip, int startrun, int stoprun) {

  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.1);

  TCanvas *c1 = new TCanvas("c1","par1",700,700);
  TProfile *ppar1 = new TProfile("par1"," ",90,0,90,-100000,100000,"");

  TCanvas *c2 = new TCanvas("c2","par2",700,700);
  TProfile *ppar2 = new TProfile("par2"," ",90,0,90,-100000,100000,"");

  TCanvas *c3 = new TCanvas("c3","par3",700,700);
  TProfile *ppar3 = new TProfile("par3"," ",90,0,90,-100000,100000,"");

  Double_t tilt_max = 0;

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

    source->cd("MyEUTelAnalysisCMSPixel");
    Double_t tilt = gettilt(inputdir,*run,chip);
    std::pair<double,double> skwpar = fitskwlin("cmsdy0vsskw");
    //std::vector<double> skwpar_pol = fitskwpol("cmsdy0vsskw",-0.1,0.1);

    ppar1->Fill(tilt,-1*skwpar.second,1);

    //ppar1->Fill(tilt,-1*skwpar_pol.at(1),1);
    //ppar2->Fill(tilt,-1*skwpar_pol.at(2),1);
    //ppar3->Fill(tilt,-1*skwpar_pol.at(3),1);

    cout << tilt << " " << skwpar.first << " " << -1*skwpar.second << endl;
    //cout << tilt << " " << -1*skwpar_pol.at(0) << " " << -1*skwpar_pol.at(1) << " " << -1*skwpar_pol.at(2) << " " << -1*skwpar_pol.at(3) << endl;
    if(tilt > tilt_max) tilt_max = tilt;
    delete source;
  }

  c1->cd();
  setStyle(ppar1,"data");
  ppar1->SetTitle(";tilt angle [#circ];skew slope");
  ppar1->GetYaxis()->SetTitleOffset(1.2);
  ppar1->GetXaxis()->SetRangeUser(0,tilt_max);
  ppar1->Draw();
  if(chip == 506) DrawFreeCMSLabels("Simulation, 150 #mum pitch",0);
  else DrawFreeCMSLabels("Simulation, 100 #mum pitch",0);

  c2->cd();
  setStyle(ppar2,"data");
  ppar2->SetTitle(";tilt angle [#circ];skew par2");
  ppar2->GetYaxis()->SetTitleOffset(1.2);
  ppar2->GetXaxis()->SetRangeUser(0,tilt_max);
  ppar2->Draw();
  if(chip == 506) DrawFreeCMSLabels("Simulation, 150 #mum pitch",0);
  else DrawFreeCMSLabels("Simulation, 100 #mum pitch",0);

  c3->cd();
  setStyle(ppar3,"data");
  ppar3->SetTitle(";tilt angle [#circ];skew par3");
  ppar3->GetYaxis()->SetTitleOffset(1.2);
  ppar3->GetXaxis()->SetRangeUser(0,tilt_max);
  ppar3->Draw();
  if(chip == 506) DrawFreeCMSLabels("Simulation, 150 #mum pitch",0);
  else DrawFreeCMSLabels("Simulation, 100 #mum pitch",0);
}
