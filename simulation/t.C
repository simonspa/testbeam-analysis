
// Daniel Pitzl, Apr 2014, MAr 2015
// result plots for turn scan
// root -l evrd.root

// .x t.C

//------------------------------------------------------------------------------
bool check( char* hs )
{
  TObject* obj = gDirectory->Get(hs);

  if( obj == NULL ) {
    cout << hs << " does not exist\n";
    return 0;
  }
  
  if( !obj->InheritsFrom( "TH1" ) ) {
    cout << hs << " is not a 1-D historgram\n";
    return 0;
  }
  cout << obj->GetTitle() << endl;
  return 1;
}

//------------------------------------------------------------------------------
double meanY( char* hs ) // for profiles only
{
  TObject* obj = gDirectory->Get(hs);
  if( !obj->InheritsFrom( "TProfile" ) ) {
    cout << hs << " is not profile plot" << endl;
    return 0;
  }
  TProfile* p = (TProfile*)obj;
  double sumw = 0;
  double sumy = 0;
  for( int i = 1; i <= p->GetNbinsX(); ++i ) { // bin 0 is underflow
    double w = p->GetBinEntries(i);
    sumw += w;
    sumy += w * p->GetBinContent(i);
  }
  if( sumw > 0.1 )
    return sumy/sumw;
  else
    return 0;
}

//------------------------------------------------------------------------------
void t()
{
  double ncol = 0;
  double lan = 0;
  double edg = 0;
  double resx = 0;

  gROOT->ProcessLine( ".x ls.C" );

  {
    char * z = "hncol";
    if( check( z ) )
      ncol = hncol->GetMean();
  }

  {
    char * z = "hqvis"; // Landau
    if( check( z ) ) {
      gROOT->ProcessLine( Form( ".x fitlang.C( \"%s\" )", z ) );

      TF1 * f1 = (TF1*)hqvis->GetListOfFunctions()->First();
      lan = f1->GetParameter(0);
    }
  }

  {
    char * z = "hpxq"; // pixel charge
    if( check( z ) ) {
      gROOT->ProcessLine( Form( ".x fitedge3.C( \"%s\" )", z ) );

      TF1 * f1 = (TF1*)hpxq->GetListOfFunctions()->First();
      edg = f1->GetParameter(0);
    }
  }

  {
    char * z = "hdxcq3"; // with Gaussian track smearing: better fit, subtract later
    if( check( z ) ) {

      resx = hdxcq3->GetRMS();
      gROOT->ProcessLine( Form( ".x fitep0.C( \"%s\", %f, %f )", z, -3*resx, 3*resx ) );

      TF1 * f1 = (TF1*)hdxcq3->GetListOfFunctions()->First();
      resx = f1->GetParameter(1);
    }
  }

  cout << gFile->GetName() << endl;
  cout << "alfa " << halfa->GetMean() << endl;
  cout
    << "  sy " << setprecision(3) << resx
    << "  col " << setprecision(4) << ncol
    << "  Lan " << setprecision(4) << lan
    << "  edg " << setprecision(3) << edg
   << endl;
}
