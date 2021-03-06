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
#include <utility>

#include "constants.h"
#include "tools.C"

using namespace std;

void fitskw_simulation() {
  cout << "inputdir!" << endl;
}

void fitskw_simulation(TString inputdir,int chip) {

  for(int tilt = 0; tilt < 91; tilt++)  {
    TString fileName;
    fileName += inputdir;
    if( !fileName.EndsWith("/") ) fileName += "/";
    fileName += "pixelav-r" + ZeroPadNumber(tilt,4) + ".out.root";

    TFile *source;
    if (!gSystem->AccessPathName( fileName )) source = TFile::Open(fileName);
    else continue;

    std::pair<double,double> skwpar;
    if(chip == 506) skwpar = fitskwlin("dx0vsskwcol");
    else skwpar = fitskwlin("dy0vsskwrow");
    std::cout << tilt << "," << skwpar.first << "," << skwpar.second << "\n";

    /*
    std::vector<double> skwpar_pol;
    if(chip == 506) skwpar_pol = fitskwpol("dx0vsskwcol",-0.1,0.1);
    else skwpar_pol = fitskwpol("dy0vsskwrow",-0.1,0.1);
    std::cout << tilt << "," << skwpar_pol.at(0) << "," << skwpar_pol.at(1) << "," << skwpar_pol.at(2) << "," << skwpar_pol.at(3) << endl;
    */

    delete source;
  }
  
}
