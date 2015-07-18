#include "../macros/tools.C"

  //------------------------------------------------------------------------------
Double_t fitSkw( Double_t *x, Double_t *par )
{
  return par[0] + par[1] * x[0];
}

//------------------------------------------------------------------------------
void fitskwlin( char* hs,int tilt, double xl=-0.1, double xr=0.1 )
{
  TH1 *h = (TH1*)gDirectory->Get(hs);

  if( h == NULL ) {
    cout << hs << " does not exist\n";
    return;
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

  tanhFcn->SetNpx(500);
  tanhFcn->SetLineWidth(4);
  tanhFcn->SetLineColor(kMagenta);

  h->Fit( "tanhFcn", "R Q", "p" );// R = range from tanhFcn

  cout << tilt << "\t" << tanhFcn->GetParameter(0) << "\t\t" << tanhFcn->GetParameter(1) << endl;
}

void fitall() {
  cout << "inputdir!" << endl;
}

void fitall(TString inputdir) {

  for(int tilt = 0; tilt < 81; tilt++)  {
    TString fileName;
    fileName += inputdir;
    if( !fileName.EndsWith("/") ) fileName += "/";
    fileName += "pixelav-r51" + ZeroPadNumber(tilt,2) + ".out.root";

    TFile *source;
    if (!gSystem->AccessPathName( fileName )) source = TFile::Open(fileName);
    else continue;

    fitskwlin("dx0vsskwcol",tilt);

    delete source;
  }
  
}
