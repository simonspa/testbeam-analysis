#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGaxis.h"
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

TGraph * getfield(TString efile) {

  TGraph *field = new TGraph();
  Double_t z, efield;
  ifstream in;
  string line;
  in.open(efile);

  while(std::getline(in,line)){
    // Skip reading comments:
    if (line[0] == '#') continue;
    if(line.empty()) continue;

    istringstream s( line );
    int i = 0;
    while (s) {
      string str;
      if(!getline( s, str, ' ' )) break;
      if(i == 0) z = atof(str.c_str()); // z position
      if(i == 1) efield = atof(str.c_str()); // field
      i++;
    }
    field->SetPoint(field->GetN(),z,efield);
  }
  in.close();
  return field;
}

void plot_efield() {

  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.1);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TCanvas *c1 = new TCanvas("c1","par0",700,700);

  TGraph * field = getfield("../simproduction/eplot_300um.out");
  TGraph * field2 = getfield("../simproduction/eplot_295um.out");
  TGraph * field3 = getfield("../simproduction/eplot_290um.out");
  TGraph * field4 = getfield("../simproduction/eplot_285um.out");

  TLegend *leg = new TLegend();
  setLegendStyle(leg);

  c1->cd();
  double max, dummy;
  field->GetPoint(0,max,dummy);
  field->GetXaxis()->SetRangeUser(0,max);
  field->SetMarkerSize(0);
  field->SetLineWidth(3);
  field->SetLineColor(kRed+1);
  field->SetTitle(";z #left[#mum#right];electric field #left[V#right]");
  leg->AddEntry(field, "Sensor 300 #mum",  "l");
  field->Draw("apl");

  field2->Draw("same lp");
  field2->SetLineWidth(3);
  field2->SetMarkerSize(0);
  leg->AddEntry(field2, "Sensor 295 #mum",  "l");

  field3->Draw("same lp");
  field3->SetLineWidth(3);
  field2->SetLineColor(kGray+3);
  field3->SetMarkerSize(0);
  leg->AddEntry(field3, "Sensor 290 #mum",  "l");

  field4->Draw("same lp");
  field4->SetLineWidth(3);
  field4->SetMarkerSize(0);
  leg->AddEntry(field4, "Sensor 285 #mum",  "l");

  DrawFreeCMSLabels("Simulation, -200 V Bias",0.045);
  leg->Draw();
}
