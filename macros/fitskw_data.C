#include "tools.C"

void fitskw_data() {
  cout << "inputdir!" << endl;
}

void fitskw_data(TString inputdir, int chip) {

  // Get all runs for given chip:
  std::vector<int> runs = getruns(inputdir,chip);
  for(std::vector<int>::iterator run = runs.begin(); run != runs.end(); run++) {

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

    cout << tilt << " " << skwpar.first << " " << skwpar.second << end;
    delete source;
  }
  
}