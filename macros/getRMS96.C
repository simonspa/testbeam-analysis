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
#include "TString.h"
#include "TROOT.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

// Returns the RMS including 96% of the histogram entries, cutting the tails:
Double_t getRMS96(char* hs) {

  bool debug = false;

  TH1 *h = (TH1*)gDirectory->Get(hs);
  if( h == NULL ){ cout << hs << " does not exist\n"; return 0; }

  // Total entries:
  double integral = h->GetEntries();
  int maxbin = h->GetMaximumBin();
  if(debug) cout << "entries=" << integral << " maxbin=" << maxbin << endl;

  double subrange_integral = h->GetBinContent(maxbin);
  int bin = 0;
  while(subrange_integral < 0.96*integral) {
    bin++;
    // Add one bin to the left:
    subrange_integral += h->GetBinContent(maxbin-bin);
    // Add one bin to the right:
    subrange_integral += h->GetBinContent(maxbin+bin);
    if(debug) cout << "subrange " << (maxbin-bin) << "-" << (maxbin+bin) << ": entries=" << subrange_integral << endl;
  }
  if(debug) cout << "subrange " << (maxbin-bin) << "-" << (maxbin+bin) << " now has " << subrange_integral << " entries, this is " << (100.0*subrange_integral)/integral << "%" << endl;

  // Correct by overshoot bin:
  subrange_integral -= h->GetBinContent(maxbin+bin);
  subrange_integral -= h->GetBinContent(maxbin-bin);
  bin--;

  int binlow = maxbin-bin;
  int binhigh = maxbin+bin;
  if(debug) cout << "subrange " << (maxbin-bin) << "-" << (maxbin+bin) << " now has " << subrange_integral << " entries, this is " << (100.0*subrange_integral)/integral << "%" << endl;

  h->GetXaxis()->SetRange(binlow,binhigh); //to restrict range to bins binlow to binhigh
  double rms96 = h->GetRMS(); //will return the RMS within the axis range

  return rms96;
}
