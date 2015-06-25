#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include "tools.C"

using namespace std;

Double_t restel = 4.8;

void gettiltruns() {
  std::cout << "Run gettiltruns(inputdir, mintilt, maxtilt=90)" << std::endl;
}

void gettiltruns(const char* inputdir, Double_t mintilt) {
  void gettiltruns(inputdir, mintilt, 90);
}

void gettiltruns(const char* inputdir, Double_t mintilt, Double_t maxtilt) {
  // Get all runs for chip 506:
  int chip = 506;

  std::stringstream allruns;
  std::vector<int> runs = getruns(inputdir,chip);
  for(std::vector<int>::iterator run = runs.begin(); run != runs.end(); run++) {
    Double_t tilt = gettilt(inputdir,*run,chip);
    if(mintilt < tilt && tilt < maxtilt) {
      cout << "run" << *run << " (chip" << chip << ") tilt " << tilt << endl;
      allruns << *run << " ";
    }
  }
  cout << allruns.str();
}
