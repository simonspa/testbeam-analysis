#include "tools.C"

void resolution_simulation() {
  cout << "inputdir!" << endl;
}

void resolution_simulation(TString inputdir, int chip) {

  bool is_rotated = false;
  if(chip == 506) is_rotated = true;

  for(int tilt = 0; tilt < 91; tilt++)  {
    TString fileName;
    fileName += inputdir;
    if( !fileName.EndsWith("/") ) fileName += "/";
    fileName += "pixelav-r00" + ZeroPadNumber(tilt,2) + ".out.root";

    TFile *source;
    if (!gSystem->AccessPathName( fileName )) source = TFile::Open(fileName);
    else continue;

    // Y Resolution in fiducial volume & Landau peak:
    Double_t res;
    if(is_rotated) res = fitep0sigma("hdx0cq3");
    else res = fitep0sigma("hdy0cq3");

    // Skew corrected:
    Double_t res_corr;
    if(is_rotated) res_corr = fitep0sigma("hdxcq3");
    else res_corr = fitep0sigma("hdycq3");

    // If the corrected is larger then the uncorrected resokution, take the RMS instead:
    // reason: after correction the desitribution looks very different from a Gauss for small angles

    //if(res_corr > res) {
    if(tilt < 15) {
      TH1 *h;
      if(is_rotated) gDirectory->GetObject("hdxcq3",h);
      else gDirectory->GetObject("hdycq3",h);
      res_corr = h->GetRMS();
      cout << "(RMS for tilt " << tilt << ")" << endl;
    }

    Double_t lanpk = fitlang("h031",18,499);

    TH1 *nc;
    if(is_rotated) gDirectory->GetObject("hcol",nc);
    else gDirectory->GetObject("hrow",nc);

    Double_t ncol = nc->GetMean();
    int nevents = nc->GetEntries();

    // FIXME what's that?
    Double_t edge = 0;
    
    std::cout << ZeroPadNumber(tilt,4) 
	      << "  " << tilt 
	      << "  " << "0.1"
	      << "  " << nevents
	      << "  " << res
	      << "  " << res_corr
	      << "  " << ncol
	      << "  " << lanpk
	      << "  " << edge
	      << std::endl;
    delete source;
  }
  
}
