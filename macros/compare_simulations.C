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

  TCanvas *c2 = new TCanvas("c2","resolution",600,600);

  gStyle->SetOptStat(0);

  std::vector<double> vtilt_gap30 = getsimulation("tilt", chip, 294, 200);
  std::vector<double> vres_gap30 = getsimulation("res", chip, 294, 200);

  std::vector<double> vtilt_dot1 = getsimulation("tilt", chip, 294, 200, true);
  std::vector<double> vres_dot1 = getsimulation("res", chip, 294, 200, true);

  TLegend *leg = new TLegend();
  setLegendStyle(leg);

  if(!vtilt_gap30.empty()) {
    c2->cd();
    TGraph *si = new TGraph( vtilt_gap30.size(), &(vtilt_gap30[0]), &(vres_gap30[0]) ); // sim
    si->SetLineColor(2);
    si->SetLineWidth(3);
    si->SetMarkerColor(2);
    si->Draw("APL");
    setStyleAndFillLegend(si,"gap30",leg);

    TGraph *si_dot1 = new TGraph( vtilt_dot1.size(), &(vtilt_dot1[0]), &(vres_dot1[0]) ); // sim
    si_dot1->SetLineColor(2);
    si_dot1->SetLineWidth(3);
    si_dot1->SetMarkerColor(2);
    si_dot1->Draw("PL"); // without axis option: overlay
    setStyleAndFillLegend(si_dot1,"dot1",leg);

    leg->Draw();
  }

}
