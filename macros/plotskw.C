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

bool cmslogo = false;

void plotskw() {
  std::cout << "Run plotskw(file)" << std::endl;
}

void plotskw(TString myfile, int chip) {

  // Set the histogram styles:
  setHHStyle(*gStyle);
  gStyle->SetTitleYOffset(1.1);

  TCanvas *c1 = new TCanvas("c1","par0",700,700);
  TProfile *ppar0 = new TProfile("par0"," ",90,0,90,-10,10,"");

  TCanvas *c2 = new TCanvas("c2","par1",700,700);
  TProfile *ppar1 = new TProfile("par1"," ",90,0,90,-2000,1000,"");

  Double_t tilt, par0, par1;

  ifstream in;
  string line;
  in.open(myfile);

  // Skip first line:
  getline(in,line);

  while(std::getline(in,line)){
    // Skip reading comments:
    if (line[0] == '#') continue;
    if(line.empty()) continue;

    istringstream s( line );
    int i = 0;
    while (s) {
      string str;
      if(!getline( s, str, ',' )) break;
      if(i == 0) tilt = atof(str.c_str()); // tilt
      if(i == 1) par0 = atof(str.c_str());
      if(i == 2) par1 = atof(str.c_str());
      i++;
    }
    cout << "tilt " << tilt << " p0 " << par0 << " p1 " << par1 << endl;
    ppar0->Fill(tilt,par0,1);
    ppar1->Fill(tilt,par1,1);
  }


  c1->cd();
  ppar0->Draw();

  c2->cd();
  setStyle(ppar1,"sim");
  if(chip == 506) ppar1->SetTitle(";tilt angle #omega [#circ];skew slope");
  else ppar1->SetTitle(";tilt angle #alpha [#circ];skew slope");
  ppar1->GetYaxis()->SetTitleOffset(1.2);
  ppar1->GetXaxis()->SetRangeUser(0,tilt);
  ppar1->Draw();
  if(chip == 506) DrawFreeCMSLabels("Simulation, 150 #mum pitch",0);
  else DrawFreeCMSLabels("Simulation, 100 #mum pitch",0);

  in.close();
}
