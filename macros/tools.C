#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>


std::vector<int> getruns(const char * inputdir, int chip) {

  std::vector<int> runs;

  void* dirp = gSystem->OpenDirectory(inputdir);
  ifstream in;
  string line;
  in.open("runlist.csv");

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

Double_t gettilt(const char * inputdir, int run, int chip) {

  Double_t tilt = -999;
  bool found_run = false;
  bool found_chip = false;

  void* dirp = gSystem->OpenDirectory(inputdir);
  ifstream in;
  string line;
  in.open("runlist.csv");

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
      if(i == 0 && run == atoi(str.c_str())) found_run = true; // Read run number
      if(i == 4 && chip == atoi(str.c_str())) found_chip = true; // Read DUT chip id
      if(i == 13 && found_run && found_chip) { tilt = atoi(str.c_str()); break; } // Store tilt value
      i++;
    }
    if(found_run && found_chip) break;
  }

  in.close();
  return tilt;
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

Double_t fitlang( char* hs ) {

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
  int ib0 = h->FindBin(18);
  int ib9 = h->FindBin(40);

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
}

Double_t fitep0sigma( char* hs ) {

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