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

#include "plotter.C"
#include "tools.C"

using namespace std;

void compare_simulations() {
  std::cout << "Run compare_simulations(chip)" << std::endl;
}

void compare_simulations(int chip) {

  TCanvas *c2 = new TCanvas("c2","resolution",700,700);

  gStyle->SetOptStat(0);
  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.1);

  int thickness = 294;
  if(chip == 506) thickness = 308;

  int threshold = 170;
  if(chip == 506) threshold = 200;

  std::vector<double> vtilt_gap30 = getsimulation("tilt", chip, thickness, threshold);
  std::vector<double> vres_gap30 = getsimulation("resskw", chip, thickness, threshold);

  std::vector<double> vtilt_dot1 = getsimulation("tilt", chip, thickness, threshold,true);
  std::vector<double> vres_dot1 = getsimulation("resskw", chip, thickness, threshold,true);

  TLegend *leg = new TLegend();
  setLegendStyle(leg);

  if(!vtilt_gap30.empty()) {
    c2->cd();
    TGraph *si = new TGraph( vtilt_gap30.size(), &(vtilt_gap30[0]), &(vres_gap30[0]) ); // sim
    setStyle(si,"sim");
    si->SetLineColor(2);
    si->SetLineWidth(3);
    si->SetMarkerColor(2);
    if(chip == 506) si->SetTitle(";tilt angle #omega [#circ];resolution y #left[#mum#right]");
    else si->SetTitle(";tilt angle #alpha [#circ];resolution x #left[#mum#right]");
    si->Draw("APL");
    setStyle(si,"dot1");
    leg->AddEntry(si, "Linear fit",  "l");

    TGraph *si_dot1 = new TGraph( vtilt_dot1.size(), &(vtilt_dot1[0]), &(vres_dot1[0]) ); // sim
    si_dot1->SetLineColor(2);
    si_dot1->SetLineWidth(3);
    si_dot1->SetMarkerColor(2);
    si_dot1->Draw("PL"); // without axis option: overlay
    setStyle(si_dot1,"gap30");
    leg->AddEntry(si_dot1, "3rd order polynomial fit",  "l");

    leg->Draw();
  }

}
