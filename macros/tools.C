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

#include "constants.h"

using namespace std;

std::vector<int> getruns(const char * inputdir, int chip, string list = "") {

  std::vector<int> runs;

  gSystem->OpenDirectory(inputdir);
  ifstream in;
  string line;
  TString name;
  name.Form("runlist-%i%s.csv",chip,list.c_str());
  std::cout << "Opening " << name << endl;
  in.open(name);

  // Skip first line:
  getline(in,line);

  while(std::getline(in,line)){
    // Skip reading comments:
    if (line[0] == '#') continue;
    if(line.empty()) continue;

    istringstream s( line );
    int i = 0;
    int run;
    while (s) {
      string str;
      if(!getline( s, str, ',' )) break;
      if(i == 0) run = atoi(str.c_str()); // Read run number
      if(i == 4 && chip == atoi(str.c_str())) { runs.push_back(run); } // Read DUT chip id
      i++;
    }
  }

  in.close();
  return runs;
}

Double_t gettilt(const char * inputdir, int run, int chip, string list = "") {

  Double_t tilt = -999;
  bool found_run = false;
  bool found_chip = false;

  gSystem->OpenDirectory(inputdir);
  ifstream in;
  string line;
  TString name;
  name.Form("runlist-%i%s.csv",chip,list.c_str());
  in.open(name);

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
      if(i == 0 && run == atoi(str.c_str())) found_run = true; // Read run number
      if(i == 4 && chip == atoi(str.c_str())) found_chip = true; // Read DUT chip id
      if(i == 14 && found_run && found_chip) { tilt = atof(str.c_str()); break; } // Store tilt value
      i++;
    }
    if(found_run && found_chip) break;
  }

  in.close();
  return tilt;
}

Double_t getalignmentrun(const char * inputdir, int run, int chip, string list = "") {

  Int_t alignrun = 0;
  bool found_run = false;
  bool found_chip = false;

  gSystem->OpenDirectory(inputdir);
  ifstream in;
  string line;
  TString name;
  name.Form("runlist-%i%s.csv",chip,list.c_str());
  in.open(name);

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
      if(i == 0 && run == atoi(str.c_str())) found_run = true; // Read run number
      if(i == 1) alignrun = atoi(str.c_str()); // Alignment run number
      if(i == 4 && chip == atoi(str.c_str()) && found_run) { found_chip = true; break; }
      i++;
    }
    if(found_run && found_chip) break;
  }

  in.close();
  return alignrun;
}

std::vector<Double_t> getalignment(const char * inputdir, int run, int chip, string list = "") {

  std::vector<Double_t> alignment;
  bool found_run = false;
  bool found_chip = false;

  gSystem->OpenDirectory(inputdir);
  ifstream in;
  string line;
  TString name;
  name.Form("runlist-%i%s.csv",chip,list.c_str());
  in.open(name);

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
      if(i == 0 && run == atoi(str.c_str())) found_run = true; // Read run number
      if(i == 4 && chip == atoi(str.c_str()) && found_run) { found_chip = true; }
      if(found_run && found_chip) {
	if(i == 7) { alignment.push_back(atof(str.c_str())); } // Gain conversion factor
	if(i == 11) { alignment.push_back(atof(str.c_str())); } // align x
	if(i == 12) { alignment.push_back(atof(str.c_str())); } // align y
	if(i == 13) { alignment.push_back(atof(str.c_str())); } // align z
	if(i == 14) { alignment.push_back(atof(str.c_str())); } // align tilt
	if(i == 15) { alignment.push_back(atof(str.c_str())); } // align turn
	if(i == 16) { alignment.push_back(atof(str.c_str())); } // align rot
      }
      i++;
    }
    if(found_run && found_chip) break;
  }

  in.close();
  return alignment;
}

Double_t getdz(const char * inputdir, int run, int chip, string list = "") {

  Double_t dz = -999;
  bool found_run = false;
  bool found_chip = false;

  gSystem->OpenDirectory(inputdir);
  ifstream in;
  string line;
  TString name;
  name.Form("runlist-%i%s.csv",chip,list.c_str());
  in.open(name);

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
      if(i == 0 && run == atoi(str.c_str())) found_run = true; // Read run number
      if(i == 4 && chip == atoi(str.c_str())) found_chip = true; // Read DUT chip id
      if(i == 13 && found_run && found_chip) { dz = atof(str.c_str()); break; } // Store tilt value
      i++;
    }
    if(found_run && found_chip) break;
  }

  in.close();
  return dz;
}

Double_t gettrim(const char * inputdir, int run, int chip, string list = "") {

  string trimstring;
  bool found_run = false;
  bool found_chip = false;

  gSystem->OpenDirectory(inputdir);
  ifstream in;
  string line;
  TString name;
  name.Form("runlist-%i%s.csv",chip,list.c_str());
  in.open(name);

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
      if(i == 0 && run == atoi(str.c_str())) found_run = true; // Read run number
      if(i == 4 && chip == atoi(str.c_str())) found_chip = true; // Read DUT chip id
      if(i == 5 && found_run && found_chip) { trimstring = str; break; } // Store value
      i++;
    }
    if(found_run && found_chip) break;
  }

  int trim = 0;
  string delimiter = "trim";
  trimstring = trimstring.substr(trimstring.find(delimiter) + delimiter.length(), 2);
  trim = atof(trimstring.c_str());
  in.close();
  return trim;
}

Double_t getTelRes(double dz) {

  Double_t dzDUT, sigma_x;

  ifstream in;
  string line;
  in.open("resolutions.db");

  // Skip first line:
  getline(in,line);

  while(std::getline(in,line)){
    // Skip reading comments:
    if (line[0] == '#') continue;
    if(line.empty()) continue;

    istringstream s( line );
    s >> dzDUT >> sigma_x;
    if(dzDUT >= dz) break;
  }

  in.close();
  return sigma_x;
}

std::vector<double> getsimulation(std::string name, int chip, int thickness=294, int threshold=200, bool dot1=false) {

  //----------------------------------------------------------------------------
  // read sim:

  TString file;
  if(dot1) file.Form("simulation/sim%i_%iskw_thr%i_pol.dat",thickness,chip,threshold);
  else file.Form("simulation/sim%i_%iskw_thr%i.dat",thickness,chip,threshold);

  cout << "try to open sim file " << file;
  
  ifstream SIMstream( file );
  if( !SIMstream ) {
    cout << ": failed" << endl;
    return std::vector<double>();
  }
  cout << ": succeed" << endl;

  // Read file by lines:
  double pi = 4*atan(1.0);
  double wt = 180/pi;

  string rl;
  double tilt;
  int run;
  int nev;
  double ry;
  double ry_skwcorr;
  double ncol;
  double lanpk;

  double turn;
  double edge;

  vector <double> stilt;
  vector <double> sry;
  vector <double> sryskw;
  vector <double> sncol;
  vector <double> slanpk;
  vector <double> slanpkn;

  while( SIMstream.good() && ! SIMstream.eof() ) {

    getline( SIMstream, rl ); // read one line  = event into string
    istringstream simrun( rl ); // tokenize string

    simrun >> run;
    simrun >> tilt; // [deg]
    simrun >> turn;
    simrun >> nev;
    simrun >> ry;
    simrun >> ry_skwcorr;
    simrun >> ncol;
    simrun >> lanpk;
    simrun >> edge;
    stilt.push_back(tilt);
    sry.push_back(ry);
    sryskw.push_back(ry_skwcorr);
    sncol.push_back(ncol);
    slanpk.push_back(lanpk);
    slanpkn.push_back(edge);
  } // while lines

  vector<double> btilttan;
  vector<double> beta;
  vector<double> bsy;
  vector<double> bsyskw;
  vector<double> bpath;
  vector<double> btant;

  for(size_t i = 0; i < stilt.size(); i++) {
    btilttan.push_back(TMath::Tan(stilt.at(i)/180*TMath::Pi()));
    beta.push_back(-TMath::Log(TMath::Tan(TMath::TwoPi()*(90-stilt.at(i))/(2*360))));
    bpath.push_back(1 / cos( stilt.at(i) / wt ));
    btant.push_back(tan( stilt.at(i) / wt ));
    bsy.push_back(sqrt( sry.at(i)*sry.at(i) - restel_sim*restel_sim )); // subtract telescope
    bsyskw.push_back(sqrt( sryskw.at(i)*sryskw.at(i) - restel_sim*restel_sim )); // subtract telescope
  }

  if(name == "tilt") return stilt;
  else if(name == "tilttan") return btilttan;
  else if(name == "eta") return beta;
  else if(name == "path") return bpath;
  else if(name == "res") return bsy;
  else if(name == "resskw") return bsyskw;
  else if(name == "ncol") return sncol;
  else if(name == "peak") return slanpk;
  else if(name == "peaknorm") return slanpkn;
  else return std::vector<double>();
}

std::vector<double> getthresholds(std::string name, int tilt, int thickness, int chip) {

  std::vector<double> vthr;
  std::vector<double> vlanpk;
  std::vector<double> vncol;
  std::vector<double> vres;

  // Read thresholds
  for(Int_t thr = 100; thr < 500; thr += 1) {

    TString myfile;
    myfile.Form("simulation/sim%i_%iskw_thr%i.dat",thickness,chip,thr);
  
    ifstream SIMstream( myfile );
    if( !SIMstream ) { continue;}

    cout << "Fetching from " << myfile << endl;
    // Read file by lines:
    string rl;
    double tlt;
    int run;
    int nev;
    double ry;
    double ry_skwcorr;
    double ncol;
    double lanpeak;
    double turn;
    double edge;

    while( SIMstream.good() && ! SIMstream.eof() ) {

      getline( SIMstream, rl ); // read one line  = event into string
      istringstream simrun( rl ); // tokenize string

      simrun >> run;
      simrun >> tlt; // [deg]
      simrun >> turn;
      simrun >> nev;
      simrun >> ry;
      simrun >> ry_skwcorr;
      simrun >> ncol;
      simrun >> lanpeak;
      simrun >> edge;

      if(tilt < (tlt+0.5)) {
	vthr.push_back(thr*0.01);
	vlanpk.push_back(lanpeak);
	vncol.push_back(ncol);
	vres.push_back(sqrt( ry*ry - restel_sim*restel_sim )); // subtract telescope
	cout << "sim tilt " << tlt << " res " << ry << " ressub " << sqrt( ry*ry - restel_sim*restel_sim ) << " electrons " << (thr*0.01) << " ncol " << ncol << endl;
	break;
      }
    } // while lines
  }
  cout << endl;
  if(name == "thr") return vthr;
  else if(name == "lanpk") return vlanpk;
  else if(name == "ncol") return vncol;
  else if(name == "res") return vres;
  else return std::vector<double>();
}

// Gauss function:
Double_t gp0Fit( Double_t *x, Double_t *par ) {

  static int nn = 0;
  nn++;
  static double dx = 0.1;
  static double b1 = 0;
  if( nn == 1 ) b1 = x[0];
  if( nn == 2 ) dx = x[0] - b1;

  //--  Mean and width:
  double t = ( x[0] - par[0] ) / par[1];

  double pi = 3.14159265358979323846;
  double aa = dx / par[1] / sqrt(2*pi);

  double pk = par[2] * aa * exp( -0.5*t*t );

  return pk + par[3];
}

Double_t fitgp0( char* hs ) {

  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ){
    cout << hs << " does not exist\n";
    return 0;
  }

  h->SetMarkerStyle(21);
  h->SetMarkerSize(0.8);
  h->SetStats(1);
  gStyle->SetOptFit(101);

  gROOT->ForceStyle();

  double dx = h->GetBinWidth(1);
  double nmax = h->GetBinContent(h->GetMaximumBin());
  double xmax = h->GetBinCenter(h->GetMaximumBin());
  double nn = 7*nmax;

  int nb = h->GetNbinsX();
  double n1 = h->GetBinContent(1);
  double n9 = h->GetBinContent(nb);
  double bg = 0.5*(n1+n9);

  double x1 = h->GetBinCenter(1);
  double x9 = h->GetBinCenter(nb);

  // create a TF1 with the range from x1 to x9 and 4 parameters
  TF1 *gp0Fcn = new TF1( "gp0Fcn", gp0Fit, x1, x9, 4 );

  gp0Fcn->SetParName( 0, "mean" );
  gp0Fcn->SetParName( 1, "sigma" );
  gp0Fcn->SetParName( 2, "area" );
  gp0Fcn->SetParName( 3, "BG" );

  gp0Fcn->SetNpx(500);
  gp0Fcn->SetLineWidth(4);
  gp0Fcn->SetLineColor(kMagenta);
  gp0Fcn->SetLineColor(kGreen);

  // set start values for some parameters:
  gp0Fcn->SetParameter( 0, xmax ); // peak position
  gp0Fcn->SetParameter( 1, 4*dx ); // width
  gp0Fcn->SetParameter( 2, nn ); // N
  gp0Fcn->SetParameter( 3, bg );

  // N: not drawing
  // Q: quiet
  // R: use specified range
  h->Fit( "gp0Fcn", "NQR", "ep" );

  return gp0Fcn->GetParameter(1);

}

std::string ZeroPadNumber(int num, int len)
{
  std::ostringstream ss;
  ss << std::setw( len ) << std::setfill( '0' ) << num;
  return ss.str();
}

Double_t fitLandauGauss( Double_t *x, Double_t *par ) {

  static int nn=0;
  nn++;
  static double xbin = 1;
  static double b1 = 0;
  if( nn == 1 ) b1 = x[0];
  if( nn == 2 ) xbin = x[0] - b1;// bin width needed for normalization

  // Landau:
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // MP shift correction:
  double mpc = par[0] - mpshift * par[1]; //most probable value (peak pos)

  //Fit parameters:
  //par[0] = Most Probable (MP, location) parameter of Landau density
  //par[1] = Width (scale) parameter of Landau density
  //par[2] = Total area (integral -inf to inf, normalization constant)
  //par[3] = Gaussian smearing

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Range of convolution integral
  double xlow = x[0] - sc * par[3];
  double xupp = x[0] + sc * par[3];
  double step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  double sum = 0;
  double xx;
  double fland;

  for( int i = 1; i <= np/2; i++ ) {
    xx = xlow + ( i - 0.5 ) * step;
    fland = TMath::Landau( xx, mpc, par[1] ) / par[1];
    sum += fland * TMath::Gaus( x[0], xx, par[3] );

    xx = xupp - ( i - 0.5 ) * step;
    fland = TMath::Landau( xx, mpc, par[1] ) / par[1];
    sum += fland * TMath::Gaus( x[0], xx, par[3] );
  }

  return( par[2] * invsq2pi * xbin * step * sum / par[3] );
}

Double_t fitlang( char* hs, double klow = 18, double khigh = 40 ) {

  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ){
    cout << hs << " does not exist\n";
    return 0;
  }

  double aa = h->GetEntries();//normalization

  // find peak:
  int ipk = h->GetMaximumBin();
  double xpk = h->GetBinCenter(ipk);
  double sm = xpk / 9; // sigma
  double ns = sm; // noise

  // fit range:
  int ib0 = h->FindBin(klow);
  int ib9 = h->FindBin(khigh);

  double x0 = h->GetBinLowEdge(ib0);
  double x9 = h->GetBinLowEdge(ib9) + h->GetBinWidth(ib9);

  // create a TF1 with the range from x0 to x9 and 4 parameters
  TF1 *fitFcn = new TF1( "fitFcn", fitLandauGauss, x0, x9, 4 );

  fitFcn->SetParName( 0, "peak" );
  fitFcn->SetParName( 1, "sigma" );
  fitFcn->SetParName( 2, "area" );
  fitFcn->SetParName( 3, "smear" );

  fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(4);
  fitFcn->SetLineColor(kMagenta);

  // set start values:
  fitFcn->SetParameter( 0, xpk ); // peak position, defined above
  fitFcn->SetParameter( 1, sm ); // width
  fitFcn->SetParameter( 2, aa ); // area
  fitFcn->SetParameter( 3, ns ); // noise

  h->Fit("fitFcn", "NQR", "ep" );// R = range from fitFcn
  return fitFcn->GetParameter(0);
}

Double_t fitfulllang( char* hs ) {

  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ){
    cout << hs << " does not exist\n";
    return 0;
  }

  double aa = h->GetEntries();//normalization

  // find peak:
  int ipk = h->GetMaximumBin();
  double xpk = h->GetBinCenter(ipk);
  double sm = xpk / 9; // sigma
  double ns = sm; // noise

  // fit range:
  int ib0 = ipk/2;
  int ib9 = h->GetNbinsX() - 1;

  double x0 = h->GetBinLowEdge(ib0);
  double x9 = h->GetBinLowEdge(ib9) + h->GetBinWidth(ib9);

  // create a TF1 with the range from x0 to x9 and 4 parameters
  TF1 *fitFcn = new TF1( "fitFcn", fitLandauGauss, x0, x9, 4 );

  fitFcn->SetParName( 0, "peak" );
  fitFcn->SetParName( 1, "sigma" );
  fitFcn->SetParName( 2, "area" );
  fitFcn->SetParName( 3, "smear" );

  fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(4);
  fitFcn->SetLineColor(kMagenta);

  // set start values:
  fitFcn->SetParameter( 0, xpk ); // peak position, defined above
  fitFcn->SetParameter( 1, sm ); // width
  fitFcn->SetParameter( 2, aa ); // area
  fitFcn->SetParameter( 3, ns ); // noise

  h->Fit("fitFcn", "NQR", "ep" );// R = range from fitFcn
  return fitFcn->GetParameter(0);
}


// Student's t function:
Double_t ep0Fit( Double_t *x, Double_t *par ) {

  static int nn = 0;
  nn++;
  static double dx = 0.1;
  static double b1 = 0;
  if( nn == 1 ) b1 = x[0];
  if( nn == 2 ) dx = x[0] - b1;

  // Mean and width:

  double bet = par[2]; // exponent

  double pk = 0;

  if( bet > 0 ) {

    double xm  = par[0];
    double sig = par[1];

    double t = ( x[0] - xm ) / sig;

    double u = pow( fabs( t/sqrt(2.0) ), bet );

    // Normalization needs Gamma function:

    double aa = dx / sig / sqrt(8.0) * bet / TMath::Gamma( 1/bet );

    pk = par[3] * aa * exp( -u );
  }

  return pk + par[4];
}

//----------------------------------------------------------------------

int fitep0( char* hs ) {

  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ){

    cout << hs << " does not exist\n";

  }

  else{
   
    h->SetMarkerStyle(21);
    h->SetMarkerSize(0.8);
    h->SetStats(1);
    gStyle->SetOptFit(101);

    gROOT->ForceStyle();

    double dx = h->GetBinWidth(1);
    double nmax = h->GetBinContent(h->GetMaximumBin());
    double xmax = h->GetBinCenter(h->GetMaximumBin());
    double nn = 7*nmax;

    int nb = h->GetNbinsX();
    double n1 = h->GetBinContent(1);
    double n9 = h->GetBinContent(nb);
    double bg = 0.5*(n1+n9);

    double x1 = h->GetBinCenter(1);
    double x9 = h->GetBinCenter(nb);
    cout << hs << ": " << x1 << " - " << x9 << endl;

    // create a TF1 with the range from x1 to x9 and 5 parameters

    TF1 *ep0Fcn = new TF1( "ep0Fcn", ep0Fit, x1, x9, 5 );

    ep0Fcn->SetParName( 0, "mean" );
    ep0Fcn->SetParName( 1, "sigma" );
    ep0Fcn->SetParName( 2, "pow" );
    ep0Fcn->SetParName( 3, "area" );
    ep0Fcn->SetParName( 4, "BG" );

    ep0Fcn->SetNpx(500);
    ep0Fcn->SetLineWidth(4);
    ep0Fcn->SetLineColor(kMagenta);
    ep0Fcn->SetLineColor(kGreen);
   
    // set start values for some parameters:

    cout << hs << " " << dx << ", " << nn << ", " << xmax << endl;

    ep0Fcn->SetParameter( 0, xmax ); // peak position
    ep0Fcn->SetParameter( 1, 4*dx ); // width
    ep0Fcn->SetParameter( 2, 3.3 ); // pow
    ep0Fcn->SetParameter( 3, nn ); // N
    ep0Fcn->SetParameter( 4, bg );
    
    h->Fit( "ep0Fcn", "R", "ep" );
    // h->Fit("ep0Fcn","V+","ep");

    h->Draw("histepsame");  // data again on top
  }
  return 0;
}

Double_t fitep0sigma( char* hs, int binlow=-999, int binhigh=999) {

  TH1 *h = (TH1*)gDirectory->Get(hs);
  if( h == NULL ){ cout << hs << " does not exist\n"; return 0; }

  double dx = h->GetBinWidth(1);
  double nmax = h->GetBinContent(h->GetMaximumBin());
  double xmax = h->GetBinCenter(h->GetMaximumBin());
  double nn = 7*nmax;

  int nb = h->GetNbinsX();
  double n1 = h->GetBinContent(1);
  double n9 = h->GetBinContent(nb);
  double bg = 0.5*(n1+n9);

  double x1, x9;
  if(binlow < -900 && binhigh > 900) {
    x1 = h->GetBinCenter(1);
    x9 = h->GetBinCenter(nb);
  }
  else {
    x1 = binlow;
    x9 = binhigh;
  }

  // create a TF1 with the range from x1 to x9 and 5 parameters
  TF1 *ep0Fcn = new TF1( "ep0Fcn", ep0Fit, x1, x9, 5 );

  ep0Fcn->SetParName( 0, "mean" );
  ep0Fcn->SetParName( 1, "sigma" );
  ep0Fcn->SetParName( 2, "pow" );
  ep0Fcn->SetParName( 3, "area" );
  ep0Fcn->SetParName( 4, "BG" );

  // Start values for some parameters:
  ep0Fcn->SetParameter( 0, xmax ); // peak position
  ep0Fcn->SetParameter( 1, 4*dx ); // width
  ep0Fcn->SetParameter( 2, 3.3 ); // pow
  ep0Fcn->SetParameter( 3, nn ); // N
  ep0Fcn->SetParameter( 4, bg );
    
  h->Fit("ep0Fcn", "Q R", "ep" );
  TF1 *fit = h->GetFunction("ep0Fcn");
  return fit->GetParameter(1);
}


// Returns the RMS including 96% of the histogram entries, cutting the tails:
Double_t getRMS96(char* hs, double truncation=96.) {

  bool debug = false;

  TH1 *h = (TH1*)gDirectory->Get(hs);
  if( h == NULL ){ cout << hs << " does not exist\n"; return 0; }

  // Total entries:
  double integral = h->GetEntries();
  int maxbin = h->GetMaximumBin();
  if(debug) cout << "entries=" << integral << " maxbin=" << maxbin << endl;

  double subrange_integral = h->GetBinContent(maxbin);
  int bin = 0;
  while(subrange_integral < truncation/100*integral) {
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


//------------------------------------------------------------------------------
Double_t fitSkw( Double_t *x, Double_t *par )
{
  return par[0] + par[1] * x[0];
}

//------------------------------------------------------------------------------
std::pair<double,double> fitskwlin( char* hs, double xl=-0.1, double xr=0.1 ) {
  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ) {
    cout << hs << " does not exist\n";
    return std::make_pair(0.,0.);
  }

  int nb = h->GetNbinsX();
  double x1 = h->GetBinCenter(1); // first
  double x9 = h->GetBinCenter(nb); // last

  if( xl > x1 && xl < x9 ) x1 = xl; // left
  if( xr > x1 && xr < x9 ) x9 = xr; // right

  // create a TF1 with the range from x1 to x9 and 3 parameters
  TF1 *tanhFcn = new TF1( "tanhFcn", fitSkw, x1, x9, 2 );
  tanhFcn->SetParName( 0, "dyoff" );
  tanhFcn->SetParName( 1, "dyslp" );

  // set start values:
  tanhFcn->SetParameter( 0, 0 ); // dy off [um]
  tanhFcn->SetParameter( 1, 99 ); // dy slope [um/skw]

  h->Fit( "tanhFcn", "R Q", "p" );// R = range from tanhFcn

  return std::make_pair(tanhFcn->GetParameter(0),tanhFcn->GetParameter(1));
}

Double_t fitSkwPol( Double_t *x, Double_t *par )
{
  return par[0] + par[1] * x[0] + par[2] * x[0]*x[0] + par[3] * x[0]*x[0]*x[0];
}

std::vector<double> fitskwpol( char* hs, double xl=-0.1, double xr=0.1) {
  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ) {
    cout << hs << " does not exist\n";
    return std::vector<double>();
  }

  int nb = h->GetNbinsX();
  double x1 = h->GetBinCenter(1); // first
  double x9 = h->GetBinCenter(nb); // last

  if( xl > x1 && xl < x9 ) x1 = xl; // left
  if( xr > x1 && xr < x9 ) x9 = xr; // right

  // create a TF1 with the range from x1 to x9 and 3 parameters
  TF1 *tanhFcn = new TF1( "tanhFcn", fitSkwPol, x1, x9, 4 );
  tanhFcn->SetParName( 0, "dyoff" );
  tanhFcn->SetParName( 1, "dyslp" );
  tanhFcn->SetParName( 2, "dypar" );
  tanhFcn->SetParName( 3, "dyhyp" );

  // set start values:
  tanhFcn->SetParameter( 0, 0 ); // dy off [um]
  tanhFcn->SetParameter( 1, 99 ); // dy slope [um/skw]
  tanhFcn->SetParameter( 2, 0 ); // dy parabola
  tanhFcn->SetParameter( 3, 0 ); // dy hyperbola


  h->Fit( "tanhFcn", "R Q", "p" );// R = range from tanhFcn

  std::vector<double> result;
  for(size_t i = 0; i < 4; i++) {
    result.push_back(tanhFcn->GetParameter(i));
  }
  return result;
}

Double_t tp0Fit( Double_t *x, Double_t *par5 ) {

  static int nn = 0;
  nn++;
  static double dx = 0.1;
  static double b1 = 0;
  if( nn == 1 ) b1 = x[0];
  if( nn == 2 ) dx = x[0] - b1;

  // Mean and width:

  double xm = par5[0];
  double t = ( x[0] - xm ) / par5[1];
  double tt = t*t;

  // exponent:

  double rn = par5[2];
  double xn = 0.5 * ( rn + 1.0 );

  // Normalization needs Gamma function:

  double pk = 0.0;

  if( rn > 0.0 ) {

    double pi = 3.14159265358979323846;
    double aa = dx / par5[1] / sqrt(rn*pi) * TMath::Gamma(xn) / TMath::Gamma(0.5*rn);

    pk = par5[3] * aa * exp( -xn * log( 1.0 + tt/rn ) );
  }

  return pk + par5[4];
}

//----------------------------------------------------------------------

int fittp0( char* hs ) {

  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ){

    cout << hs << " does not exist\n";

  }

  else{
   
    h->SetMarkerStyle(21);
    h->SetMarkerSize(0.8);
    h->SetStats(1);
    gStyle->SetOptFit(101);

    gROOT->ForceStyle();

    double dx = h->GetBinWidth(1);
    double nmax = h->GetBinContent(h->GetMaximumBin());
    double xmax = h->GetBinCenter(h->GetMaximumBin());
    double nn = 7*nmax;

    int nb = h->GetNbinsX();
    double n1 = h->GetBinContent(1);
    double n9 = h->GetBinContent(nb);
    double bg = 0.5*(n1+n9);

    double x1 = h->GetBinCenter(1);
    double x9 = h->GetBinCenter(nb);
    cout << hs << ": " << x1 << " - " << x9 << endl;

    // create a TF1 with the range from x1 to x9 and 5 parameters

    TF1 *tp0Fcn = new TF1( "tp0Fcn", tp0Fit, x1, x9, 5 );

    tp0Fcn->SetParName( 0, "mean" );
    tp0Fcn->SetParName( 1, "sigma" );
    tp0Fcn->SetParName( 2, "nu" );
    tp0Fcn->SetParName( 3, "area" );
    tp0Fcn->SetParName( 4, "BG" );

    tp0Fcn->SetNpx(500);
    tp0Fcn->SetLineWidth(4);
    tp0Fcn->SetLineColor(kMagenta);
    tp0Fcn->SetLineColor(kGreen);
   
    // set start values for some parameters:

    cout << hs << " " << dx << ", " << nn << ", " << xmax << endl;

    tp0Fcn->SetParameter( 0, xmax ); // peak position
    tp0Fcn->SetParameter( 1, 4*dx ); // width
    tp0Fcn->SetParameter( 2, 2.2 ); // nu
    tp0Fcn->SetParameter( 3, nn ); // N
    tp0Fcn->SetParameter( 4, bg );
    
    h->Fit( "tp0Fcn", "R", "ep" );
    // h->Fit("tp0Fcn","V+","ep");

    h->Draw("histepsame");  // data again on top
  }
}

double fittp0sigma( char* hs ) {

  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ){ cout << hs << " does not exist\n"; return 0; }

  double dx = h->GetBinWidth(1);
  double nmax = h->GetBinContent(h->GetMaximumBin());
  double xmax = h->GetBinCenter(h->GetMaximumBin());
  double nn = 7*nmax;

  int nb = h->GetNbinsX();
  double n1 = h->GetBinContent(1);
  double n9 = h->GetBinContent(nb);
  double bg = 0.5*(n1+n9);

  double x1 = h->GetBinCenter(1);
  double x9 = h->GetBinCenter(nb);

  // create a TF1 with the range from x1 to x9 and 5 parameters
  TF1 *tp0Fcn = new TF1( "tp0Fcn", tp0Fit, x1, x9, 5 );

  tp0Fcn->SetParName( 0, "mean" );
  tp0Fcn->SetParName( 1, "sigma" );
  tp0Fcn->SetParName( 2, "nu" );
  tp0Fcn->SetParName( 3, "area" );
  tp0Fcn->SetParName( 4, "BG" );
   
  // set start values for some parameters:
  tp0Fcn->SetParameter( 0, xmax ); // peak position
  tp0Fcn->SetParameter( 1, 4*dx ); // width
  tp0Fcn->SetParameter( 2, 2.2 ); // nu
  tp0Fcn->SetParameter( 3, nn ); // N
  tp0Fcn->SetParameter( 4, bg );
    
  h->Fit( "tp0Fcn", "Q R", "ep" );
  // h->Fit("tp0Fcn","V+","ep");
  TF1 *fit = h->GetFunction("tp0Fcn");
  return fit->GetParameter(1);
}
