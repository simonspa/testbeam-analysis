
// Daniel Pitzl, 2013, Apr 2014
// read pixelav events, add thr, noi, smr, make plots

// evrd -t 140 -c 205 pixelav-r0094.out  tilt 19.5, thick 300

// evrd -t 110 -c 205 pixelav-r0095.out  tilt 19.5, thick 285, for Tsunami corrected data

// evrd -t 120 -c 205 pixelav-r0106.out  tilt 20.0, thick 285, for 6692

// evrd -t 120 -c 205 294-19p5.out  tilt 19.5 thick 294
// evrd -t 160 -c 500 294-19p5.out  tilt 19.5 thick 294

// evrd -t 120 -c 205 pixelav-r0114.out // 0 deg, gainls

// evrd -t 215 -c 202 pixelav-r0092.out
// evrd -t 215 -c 203 pixelav-r0093.out

// evrd -t 180 pixelav-r0173.out

//evrd -c 1203 -t 180 pixelav-r0201.out
//evrd -c 1203 -t 180 pixelav-r1001.out

//evrd -c 506 -t 150 pixelav-r5027.out

#include <stdlib.h> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <fstream> // files
#include <vector>
#include <cmath>
#include <math.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom1.h>

using namespace std;

#include "pixelForReadout.h" //struct pixel, struct cluster

pixel pb[4160]; // global declaration: vector of pixels with hit
int fNHit; // global
int fCluCut = 1; // clustering: 1 = no gap (15.7.2012)

// ----------------------------------------------------------------------
vector<cluster> getClus()
{
  // returns clusters with local coordinates
  // decodePixels should have been called before to fill pixel buffer pb 
  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  vector<cluster> v;
  if( fNHit == 0 ) return v;

  int* gone = new int[fNHit];

  for( int i = 0; i < fNHit; i++ )
    gone[i] = 0;

  int seed = 0;

  while( seed < fNHit ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( int i = 0; i < fNHit; i++ ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); p++ ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut) 
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list of clusters:

    c.sumA = 0;
    c.charge = 0;
    c.size = 0;
    c.col = 0;
    c.row = 0;
    double sumQ = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  p++ ) {
      c.sumA += p->ana; // Aout
      double Qpix = p->anaVcal; // calibrated [ke]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      c.charge += Qpix;
      //if( Qpix > 20 ) Qpix = 20; // DP 25.8.2013, cut tail. tiny improv
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
    }

    c.size = c.vpix.size();

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( fabs( c.charge ) > 0.1 ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with zero charge" << endl;
    }

    v.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( (++seed < fNHit) && gone[seed] );

  } // while over seeds

  // nothing left,  return clusters

  delete gone;
  return v;
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give file name" << endl;
    return 1;
  }

  // file name = last argument:

  string evFileName( argv[argc-1] );

  cout << "try to open  " << evFileName;

  ifstream evFile( argv[argc-1] );

  if( !evFile ) {
    cout << " : failed " << endl;
    return 2;
  }

  cout << " : succeed " << endl;

  // extract run number from pixelav-r1001.out

  string::size_type idx = evFileName.find( ".out" );

  if( idx == string::npos ) {
    cout << evFileName << " does not end wih .out. end" << endl;
    return 2;
  }

  string runName = evFileName.substr( idx-4, idx ); // 1001
  int run = atoi( runName.c_str() ); // c++11: stoi( runName )

  cout << "run " << runName << " = " << run << endl;

  // further arguments:

  double thr = 150; // [dekae]
  int chip = 205;

  double cut0 = 0.0;
  double cut9 = 0.0;

  for( int i = 1; i < argc; i++ ) {

    if( !strcmp( argv[i], "-t" ) )
      thr = atoi( argv[++i] ); // dekae

    if( !strcmp( argv[i], "-c" ) )
      chip = atoi( argv[++i] );

    if( !strcmp( argv[i], "-a" ) )
      cut0 = atof( argv[++i] );

    if( !strcmp( argv[i], "-e" ) )
      cut9 = atof( argv[++i] );

  } // argc

  cout << "thr = " << thr*10 << " e" << endl;

  //double thrke = 0.01 * thr; // [ke]

  // electronics noise:  (** default Dec 2013 for chip 205

  double noi = 18; //** [dekae], data noise, OK (chip 205, r0440)

  // cross talk: makes clusters larger

  double fct = 0.00; // dig after tsunami correction
  // good overall description, need thickness 300 um
  //double fct = 0.025; // run  95 for chip 205 run 6692
  // cmsnrowvsym not well described: data dip to 1.45, sim dip to 1.7
  //double fct = 0.04; // run 95: thick 285, thr 160, large clusters too large

  // gain factor:

  double gf = 1.0;
  if( chip == 202 ) gf = 1.02; // fit for run 80 (thick 300, tilt 19.5)
  if( chip == 202 ) gf = 0.98; // fit for run 86 (thick 285, tilt 22.0)

  if( chip == 1203 ) gf = 0.92; // 203i at zero fluence corner

  //if( chip == 205 ) gf = 1.01; // fit for run 440 (thick 300, tilt 19.5)
  if( chip == 205 ) gf = 1.04; // fit for run 95 (thick 285, tilt 19.5)
  //if( chip == 205 ) gf = 1.05; // fit for run 294-14p5

  if( chip ==-205 ) gf = 1.03; // fit for run 114 (thick 285, tilt  1.5)

  if( chip == 500 ) gf = 1.04; // 294-19p5
  if( chip == 506 ) gf = 1.04; // 294-19p5

  gf *= 1-fct; // keep normalization after crosstalk

  // multiplicative gain smearing for Landau peak width:

  //double gsmr = 0.10; // Landau peak too broad
  //double gsmr = 0.08; // Landau peak OK, tail too high
  double gsmr = 0.06; // Landau OK
  //double gsmr = 0.00; // Landau peak too narrow

  // threshold smearing:

  //double tsmr = 0; // [dekae] 
  double tsmr = 10; // [dekae] measured
  //double tsmr = 20; // [dekae] npx OK

  // ADC smearing for pixel charge left edge

  //double asmr = 40; // [dekae] pxq for 205 run 10928
  double asmr = 50; // [dekae] pxq for 205 run 6692

  // telescope track smearing:

  //double smr = 9.3; // [um] DUTz 117 mm at 4.4 GeV (turn)
  //double smr = 4.1; // [um] chip 202 ex105-32 at 6 GeV, rmsdy 6.5 (chip 205, r0440)
  //double smr = 6.0; // ncolvsxm too broad, rmsdy 7.1 (chip 205, r0440)
  //double smr = 5.0; // ncolvsxm a little too narrow (chip 205, r0440)
  double smr = 5.5; // ncolvsxm OK, rmsdy 7.0 (chip 205, r0440)

  if( chip == 202 ) smr = 5.4; // nominal for runs 4732-4740
  if( chip == 202 ) smr = 6.6; // sim 93 for runs 4732-4740

  if( chip ==  203 ) smr = 5.4; // nominal for runs 4964-4980
  if( chip ==  203 ) smr = 7.0; // sim 92 for runs 4964-4980
  if( chip == 1203 ) smr = 4.2; // 203i 19 deg run 11191-
  if( chip == 1203 ) smr = 5.0; // 203i 72 deg run 11251-

  if( chip == 205 ) smr = 4.1; // chip205 nominal run 10891
  if( chip == 205 ) smr = 5.0; // chip205, r95
  if( chip == 205 ) smr = 4.8; // chip205, r101
  if( chip == -205 ) smr = 4.1; // run 12218

  if( chip == 500 ) smr = 5.9; // Jan 2015

  if( chip == 506 ) smr = 5.6; // Feb 2015

  double wbg = 0.02; // flat background for resolution profiles
  //double wbg = 0.00; // test: no BG

  double pi = 4 * atan(1);
  double wt = 180/pi;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TString histfilename;
  histfilename.Form("%s.root",evFileName.c_str());
  TFile* histoFile = new TFile( histfilename, "RECREATE" );

  // book histos:

  TH1D hxmid( "hxmid", "x mid;x [um];events", 150, -75,  75 );
  TH1D hxmod( "hxmod", "x mid;x [um];events", 150, -75,  75 );
  TH1D humod( "humod", "x mid;x [um];events", 150, -75,  75 );
  TH1D hymid( "hymid", "y mid;y [um];events", 100, -50,  50 );
  TH1D hzmid( "hzmid", "z mid;z [um];events", 285,   0, 285);

  TH1D hux( "hux", "x direction;ux;events", 100, 0, 1 );
  TH1D huy( "huy", "y direction;uy;events", 100, 0, 1 );
  TH1D huz( "huz", "z direction;uz;events", 100, 0, 1 );

  TH1D halfa( "halfa", "turn angle;turn angle #alpha [^{o}];events", 90, 0, 90 );
  TH1D hbeta( "hbeta", "tilt angle;tilt angle #beta [^{o}];events", 90, 0, 90 );

  TH1D hpxqa( "hpxqa", "raw pixel charge;pixel charge [ke];pixels", 100, 0, 10 );
  TH1D hpxqt( "hpxqt", "raw pixel charge;pixel charge [ke];pixels", 100, 0, 10 );
  TH1D hpxqs( "hpxqs", "raw pixel charge;pixel charge [ke];pixels", 100, 0, 10 );

  TH1D hq( "hq", "generated e-h pairs;generated e-h pairs [ke];events", 150, 0, 150 );
  TH1D hqvis( "hqvis", "cluster charge;cluster charge [ke];clusters", 150, 0, 150 );
  TH1D hqvis1( "hqvis1", "cluster charge no bias dot;cluster charge [ke];clusters away from bias dot", 150, 0, 150 );
  TH1D hqnrm( "hqnrm", "normal cluster charge;normal cluster charge [ke];clusters", 150, 0, 150 );
  TH1D hnall( "hnall", "pixels with charge;generated pixels;clusters", 21, -0.5, 20.5 );
  TH1D hnvis( "hnvis", "pixels above threshold;pixels above threshold;clusters", 21, -0.5, 20.5 );
  TH1D hncol( "hncol", "cols per cluster;cluster size [cols];clusters", 21, -0.5, 20.5 );
  TH1D hnrow( "hnrow", "rows per cluster;cluster size [rows];clusters", 21, -0.5, 20.5 );

  TProfile pxqvsd( "pxqvsd", "pixel q vs track depth;track depth [mm];<pixel q> [ke]", 200, -0.2, 0.2, 0, 99 );
  TProfile weibull( "weibull", "Weibull;q [ke];PH [ADC]", 100, 0, 100, -100, 999 );

  TProfile nrowvsy( "nrowvsy", "rows/cluster;y mid [um];<row/cluster>", 200, 0, 200, 0, 11 );

  TProfile ncolvsxmid( "ncolvsxmid", "cols/cluster;x track [um];<col/cluster>", 150, 0, 300, 0, 11 );
  TProfile ncolvsxm( "ncolvsxm", "cols/cluster;x track [um];<col/cluster>", 150, 0, 300, 0, 11 );
  TProfile nrowvsxm( "nrowvsxm", "rows/cluster;x track [um];<row/cluster>", 150, 0, 300, 0, 11 );
  TProfile nrowvsym( "nrowvsym", "rows/cluster;y track [um];<row/cluster>", 100, 0, 200, 0, 11 );
  TProfile2D npxvsxmym( "npxvsxmym", "pixels/cluster;x track [um];y track [um];<pixels/cluster>",
		  60, 0, 300, 40, 0, 200, 0, 9 );

  TProfile qvsxm( "qvsxm", "charge/cluster;x mid [um];<charge/cluster> [ke]", 150, 0, 300, 0, 99 );
  TProfile qvsym( "qvsym", "charge/cluster;y mid [um];<charge/cluster> [ke]", 100, 0, 200, 0, 99 );
  TProfile2D qvsxmym( "qvsxmym", "charge/cluster;x mid [um];y mid [um];<charge/cluster> [ke]",
		  60, 0, 300, 40, 0, 200, 0, 99 );

  TH1D h003( "h003", "clusters per event;clusters/event;triggers", 10,  0.5, 10.5 );
  TH1D h030( "h030", "cluster size;pixels/cluster;cluster", 21, -0.5, 20.5 );
  TH1D h031( "h031", "cluster charge;cluster charge [ke];clusters", 200, 0, 100 );

  TH1D h036( "h036", "cluster center col;cluster center [column];clusters", 21, -0.5, 20.5 );
  TH1D h037( "h037", "cluster center row;cluster center [row];clusters",    11, -0.5, 10.5 );

  TH1D hcol( "hcol", "col/clus;col/cluster;cluster", 11, -0.5, 10.5 );
  TH1D hrow( "hrow", "row/clus;row/cluster;cluster", 21, -0.5, 20.5 );
  TH1D hq1( "hq1", "2-pix cluster Q1;2-pixel cluster Q1 [ke];clusters", 200, 0, 100 );
  TH1D hq2( "hq2", "2-pix cluster Q2;2-pixel cluster Q2 [ke];clusters", 200, 0, 100 );
  TH1D heta( "heta", "eta 2-row cluster;2-row eta;cluster", 100, -1, 1 );

  TH1D hdx0( "hdx0", "dx;dx [um];clusters", 250, -250, 250 );
  TH1D hdy0( "hdy0", "dy;dy [um];clusters", 300, -100, 200 );

  TH1D hdx( "hdx", "dx;dx [um];clusters", 250, -250, 250 );
  TH1D hdy( "hdy", "dy;dy [um];clusters", 300, -100, 200 );

  TH1D skwcol( "skwcol", "skwcol;cluster col skw;clusters", 80, -0.2, 0.2 );
  TProfile skwcolvsxm( "skwcolvsxm", "skwcolvsxm;x track [um];cluster col skw", 150, 0, 80, -0.2, 0.2);
  TProfile dxvsskwcol( "dxvsskwcol", "dxvsskwcol;cluster col skw;DUT cluster - track x", 80, -0.2, 0.2, -60, 60);
  TProfile dx0vsskwcol( "dx0vsskwcol", "dx0vsskwcol;cluster col skw;DUT cluster - track x", 80, -0.2, 0.2, -60, 60);

  TH1D hdxc( "hdxc", "dx;dx [um];clusters", 250, -250, 250 );
  TH1D hdyc( "hdyc", "dy;dy [um];clusters", 300, -100, 200 );
  TH1D hdxcq0( "hdxcq0", "dx;#Deltax [um];clusters", 250, -250, 250 );
  TH1D hdycq0( "hdycq0", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdxcq1( "hdxcq1", "dx;#Deltax [um];clusters", 250, -250, 250 );
  TH1D hdycq1( "hdycq1", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdxcq2( "hdxcq2", "dx;#Deltax [um];clusters", 250, -250, 250 );
  TH1D hdycq2( "hdycq2", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdxcq3( "hdxcq3", "dx;#Deltax [um];clusters", 250, -250, 250 );
  TH1D hdx0cq3( "hdx0cq3", "dx0;#Deltax0 [um];clusters", 250, -250, 250 );
  TH1D hdycq3( "hdycq3", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdycq3d( "hdycq3d", "dytr with bg no dot;#Deltaytr [um];clusters", 300, -100, 200 );
  TH1D hdycq5( "hdycq5", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdycq450( "hdycq450", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdycq375( "hdycq375", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdycq300( "hdycq300", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdycq225( "hdycq225", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdycq187( "hdycq187", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdycq150( "hdycq150", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdycq075( "hdycq075", "dy;#Deltay [um];clusters", 300, -100, 200 );
  TH1D hdymid( "hdymid", "dy;#Deltay [um];clusters", 300, -100, 200 );

  TProfile rmsxvsxm( "rmsxvsxm", "x resolution vs x;x track [um];rms(#Deltax) [#mum]", 150, 0, 300, 0, 150 );
  TProfile rmsyvsxm( "rmsyvsxm", "y resolution vs x;x track [um];rms(#Deltay) [#mum]", 150, 0, 300, 0, 100 );
  TProfile rmsyvsym( "rmsyvsym", "y resolution vs y;y track [um];rms(#Deltay) [#mum]", 100, 0, 200, 0, 100 );
  TProfile2D rmsyvsxmym( "rmsyvsxmym", "y resolution vs x-y;x track [um];y track [um];rms(#Deltay) [#mum]", 60, 0, 300, 40, 0, 200, 0, 100 );

  TProfile dyvsxm( "dyvsxm", "y shift vs x;x track [um];shift(#Deltay) [#mum]", 150, 0, 300, -100, 100 );
  TProfile dyvsym( "dyvsym", "y shift vs y;y track [um];shift(#Deltay) [#mum]", 100, 0, 200, -100, 200 );

  TProfile dx0vsxm( "dx0vsxm", "x0 shift vs x;x track [um];shift(#Deltax) [#mum]", 150, 0, 300, -100, 100 );
  TProfile dxvsxm( "dxvsxm", "x shift vs x;x track [um];shift(#Deltax) [#mum]", 150, 0, 300, -100, 100 );
  TProfile dxvsym( "dxvsym", "x shift vs y;y track [um];shift(#Deltax) [#mum]", 100, 0, 200, -100, 200 );

  TProfile etavsym( "etavsym", "eta vs track y;y track [um];2-row <eta>", 100, 0, 200, -1, 1 );
  TProfile q0vsym( "q0vsym", "Qleft vs track y;y track [um];2-row <Qleft> [ke]", 100, 0, 200, 0, 99 );
  TProfile q9vsym( "q9vsym", "Qright vs track y;y track [um];2-row <Qright> [ke]", 100, 0, 200, 0, 99 );
  TProfile q2vsym( "q2vsym", "Q2 vs track y;y track [um];2-row <Q2> [ke]", 100, 0, 200, 0, 99 );
  TProfile rmsyvsq2( "rmsyvsq2", "y resolution vs q2;Q_{2} [ke];MAD(#Deltay) [#mum]", 100, 0, 20, 0, 100 );
  TProfile rmsyvseta( "rmsyvseta", "y resolution vs eta;eta;MAD(#Deltay) [#mum]", 100, -1, 1, 0, 100 );

  TProfile rmsxvsq( "rmsxvsq", "x resolution vs charge;cluster charge [ke];MAD(#Deltax) [#mum]", 140, 10, 150, 0, 150 );
  TProfile rmsyvsq( "rmsyvsq", "y resolution vs charge;cluster charge [ke];MAD(#Deltay) [#mum]", 150,  0, 150, 0, 100 );

  double qvec[999];
  double dq = 1.01;
  qvec[0] = 0;
  int ii = 1;
  do{
    qvec[ii] = qvec[ii-1] + pow( dq, ii );
    ii++;
  }
  while( qvec[ii-1] < 150 && ii < 999 );
  //for( int jj = 0; jj < ii; ++jj ) cout << jj << "  " << qvec[jj] << endl;

  TProfile rmsyvsqv( "rmsyvsqv", "y resolution vs charge;cluster charge [ke];MAD(#Deltay) [#mum]", ii-1, qvec, 0, 200 );
  TH1D hqnrmv( "hqnrmv", "normal cluster charge;normal cluster charge [ke];clusters", ii-1, qvec );
  TH1D hqvis1v( "hqvis1v", "cluster charge no bias dot;cluster charge [ke];clusters away from bias dot", ii-1, qvec );

  TProfile npxvsq( "npxvsq", "cluster size vs charge;cluster charge [ke];<cluster size> [pixels]", 140, 10, 150, 0, 20 );
  TProfile pxqvsq( "pxqvsq", "pixel charge vs cluster charge;cluster charge [ke];<pixel charge> [ke]", 140, 10, 150, 0, 99 );
  TProfile npxvsqv( "npxvsqv", "cluster size vs charge;cluster charge [ke];<cluster size> [pixels]", ii-1, qvec, 0, 20 );
  TProfile pxqvsqv( "pxqvsqv", "pixel charge vs cluster charge;cluster charge [ke];<pixel charge> [ke]", ii-1, qvec, 0, 99 );
  TProfile pxqvsxm( "pxqvsxm", "pixel q vs track x;x track [um];<pixel q> [ke]", 150, 0, 300, 0, 99 );
  TProfile pxqvsym( "pxqvsym", "pixel q vs track y;y track [um];<pixel q> [ke]", 100, 0, 200, 0, 99 );
  TH1D hpxq( "hpxq", "pixel charge;pixel charge [ke];pixels", 100, 0, 25 );
  TH1D hpxq1( "hpxq1", "pixel charge 1-pix;pixel charge [ke];pixels", 100, 0, 25 );
  TH1D hpxq2( "hpxq2", "pixel charge 2-pix;pixel charge [ke];pixels", 100, 0, 25 );
  TH1D hpxq3( "hpxq3", "pixel charge 3-pix;pixel charge [ke];pixels", 100, 0, 25 );

  TRandom1 * myRandom = new TRandom1(); // 1 = Ranlux

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // header line:

  string s;
  getline( evFile, s );
  cout << s << endl;

  // geometry:

  double pitchx;
  double pitchy;
  double thick;
  evFile >> pitchx;
  evFile >> pitchy;
  evFile >> thick;
  cout << "pixel volume " << pitchx << " x " << pitchy << " x " << thick << endl;

  //const int Ncol = 21; // data array
  //const int Nrow =  7; // static float pixel[2][21][7];

  // from run 1001: shallow pixelav

  const int Ncol = 13; // data array
  const int Nrow = 21; // static float pixel[2][13][21];

  if( run < 1001 && ( Ncol != 21 || Nrow != 7 ) ) {
    cout << "wrong Ncol, Nrow for run " << run << endl;
    return 3;
  }

  if( run > 1000 && ( Ncol != 13 || Nrow != 21 ) ) {
    cout << "wrong Ncol, Nrow for run " << run << endl;
    return 3;
  }

  // gain factor:

  gf *= 285 / thick; // reduce to nominal thickness

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  int iev = 0;
  bool ldb = 0;

  while( ! evFile.eof() ) {

    iev++;
    if( ldb ) cout << endl << "event " << iev << endl;
    if( !(iev%1000) ) cout << "event " << iev << endl;

    // read track position at z = zmid:

    double xmid, ymid, zmid;
    evFile >> xmid; // -75..75
    evFile >> ymid; // -50..50
    evFile >> zmid;
    if( ldb ) cout << "mid " << xmid << ", " << ymid << ", " << zmid << endl;
    hxmid.Fill( xmid );
    hymid.Fill( ymid );
    hzmid.Fill( zmid );

    // read track direction:

    double ux, uy, uz;
    evFile >> ux;
    evFile >> uy;
    evFile >> uz;
    if( ldb ) cout << "direction " << ux << ", " << uy << ", " << uz << endl;
    hux.Fill( ux );
    huy.Fill( uy );
    huz.Fill( uz );
    double alfa = atan2( ux, uz ); // turn
    double beta = atan2( uy, uz ); // tilt
    halfa.Fill( alfa*wt );
    hbeta.Fill( beta*wt );
    //double length = sqrt( 1 + ux*ux + uy*uy ); // path length / thickness
    double length = 1/cos(alfa)/cos(beta);

    // telescope track smearing:

    double rx, ry;
    myRandom->Rannor( rx, ry ); // get two normal randoms
    double xtr = xmid + smr*rx/cos(alfa);
    double ytr = ymid + smr*ry/cos(beta);

    // reduce to one pixel: (fmod does not like negatives)

    double xmod = fmod( xtr + 9.5*pitchx, pitchx ) - 0.5*pitchx;
    double umod = fmod(-xtr + 9.5*pitchx, pitchx ) - 0.5*pitchx; // flipped bias dot
    double ymod = fmod( ytr + 9.5*pitchy, pitchy ) - 0.5*pitchy;

    hxmod.Fill( xmod );
    humod.Fill( umod );

    // read number of generated e-h pairs:

    int Neh;
    evFile >> Neh;
    if( ldb ) cout << "e-h pairs " << Neh << endl;
    hq.Fill( Neh*1E-3 );

    // read pixel charge map:

    pixel all0[Ncol*Nrow];
    int Nall = 0;

    int colmin = 99;
    int colmax = -1;

    for( int jj = 0; jj < Nrow; ++jj ) // 0..20, center 10
      for( int ii = 0; ii < Ncol; ++ii ) { // 0..12, center 6

	double q;
	evFile >> q; // [10e]

	if( q > thr ) {
	  hpxqa.Fill( q*0.01 ); // [ke]
	  if( ii < colmin ) colmin = ii;
	  if( ii > colmax ) colmax = ii;
	}

	all0[Nall].ana = q; // [dekae]
	all0[Nall].col = ii;
	all0[Nall].row = jj;

	Nall++;
      }

    // row cross talk: (move after threshold?)

    pixel all1[Ncol*Nrow];

    for( int ii = 0; ii < Nall; ++ii ) {
      all1[ii].col = all0[ii].col;
      all1[ii].row = all0[ii].row;
      all1[ii].ana = all0[ii].ana;
    }

    for( int ii = 0; ii < Nall-1; ++ii )
      for( int jj = ii+1; jj < Nall; ++jj ) // all pairs
	if( all0[ii].col == all0[jj].col && // same col
	    all0[ii].row + 1 == all0[jj].row ) { // adjacent rows
	  all1[ii].ana += fct * all0[jj].ana; // mutual cross talk
	  all1[jj].ana += fct * all0[ii].ana;
	}

    int Nvis = 0;
    double qvis = 0;

    double qrow[Nrow];
    for( int row = 0; row < Nrow; ++row ) qrow[row] = 0;

    double keep1st = myRandom->Rndm(); // 0..1
    double keeplst = myRandom->Rndm(); // 0..1

    for( int ii = 0; ii < Nall; ++ii ) {

      if( all1[ii].col == colmin && keep1st < cut0 ) continue; // kill 1st col
      if( all1[ii].col == colmax && keeplst < cut9 ) continue; // kill lst col

      double q = all1[ii].ana;

      q *= gf; // gain

      // add noise: cluster size

      q += myRandom->Gaus( 0, noi ); // [dekae]

      // threshold smearing:

      double t = thr + myRandom->Gaus( 0, tsmr );

      // multiplicative gain smearing: A = g*Q
      // simulate gain saturation (Weibull)?

      // Weibull response (gainls2ps)
      // Weibull x = large Vcal = 350 e = 35 dekae

      const int npar = 3;
      double par[npar];
      par[0] =-149983.44; // horizontal offset
      par[1] = 150016.56; // width, x scaled
      par[2] = 3127.541; // power

      double xt = ( q/35 - par[0] ) / par[1];
      double lt = log(xt);
      double tp = exp( par[2]*lt ); // t^p
      double c = 1 - exp(-tp);

      weibull.Fill( q*0.01, c ); // 0.4 .. 1

      //q *= myRandom->Gaus( 1, gsmr ); // [dekae]

      q += (c-0.4)*2000/0.4 * myRandom->Gaus( 0, gsmr ); // [dekae]

      if( q > t ) { // smeared threshold

	hpxqt.Fill( q*0.01 ); // [ke]

	// ADC smearing: chip 205, run 10928

	q += myRandom->Gaus( 0, asmr ); // [dekae]

	hpxqs.Fill( q*0.01 ); // [ke]

	qvis += q;

	pb[Nvis].ana = q; // [dekae]
	pb[Nvis].col = all1[ii].col;
	pb[Nvis].row = all1[ii].row;
	pb[Nvis].anaVcal = q*0.01; // [ke]

	qrow[all1[ii].row] += q*0.01; // [ke] project onto rows in road

	Nvis++;

      } // thr

    } // all

    for( int row = 0; row < Nrow; ++row ) {

      double ypix = ( row + 0.5 - 0.5*Nrow ) * pitchy; // track enters in central pixel
      double depth = 0;
      if( beta*wt > 1 )
	depth = ( ypix - ytr ) / tan(beta); // depth [-thick/2,thick/2]

      pxqvsd.Fill( depth*1E-3, qrow[row] );

    } // rows

    double qnrm = qvis / length;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // find clusters in pb:

    fNHit = Nvis; // for cluster search

    vector<cluster> clust = getClus();	    

    h003.Fill( clust.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // look at clusters:

    for( vector<cluster>::iterator c = clust.begin(); c != clust.end(); c++ ) {

      h030.Fill( c->size );
      h031.Fill( c->charge );

      h036.Fill( c->col );
      h037.Fill( c->row );

      // pix in clus:

      int colmin = 99;
      int colmax = -1;
      int rowmin = 99;
      int rowmax = -1;

      double qcol[Ncol];
      for( int icol = 0; icol < Ncol; ++icol ) qcol[icol] = 0;

      double qrow[Nrow];
      for( int row = 0; row < Nrow; ++row ) qrow[row] = 0;

      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ ) {

	if( px->col < colmin ) colmin = px->col;
	if( px->col > colmax ) colmax = px->col;
	if( px->row < rowmin ) rowmin = px->row;
	if( px->row > rowmax ) rowmax = px->row;

	qcol[px->col] += fabs(px->anaVcal); // [ke] project cluster onto cols
	qrow[px->row] += fabs(px->anaVcal); // [ke] project cluster onto rows
	
      } // pixel loop

      int ncol = colmax - colmin + 1;
      int nrow = rowmax - rowmin + 1;

      hcol.Fill( ncol );
      hrow.Fill( nrow );

      // skew calculation (3rd moment):
      double skw = 0;
      // Sum third powers of x-x_cog:
      for(int col = colmin; col <= colmax; col++) { skw += pow((col - c->col),3)*qcol[col]; }
      // Normalize to total charge and cluster length ^3:
      skw /= (c->charge*pow(ncol,3));
      skw *= 8;

      // eta-algo in rows:

      double q1 = 0;
      double q2 = 0;
      //int col1 = 99;
      //int col2 = 99;
      int row1 = 99;
      int row2 = 99;
      double sumq = 0;
      double sumrow = 0;
      double sumrow2 = 0;

      for( int row = rowmin; row <= rowmax; ++row ) {

	double q = qrow[row]; // [ke]
	sumq += q;
	sumrow += row*q;
	sumrow2 += row*row*q;

	if( q > q1 ) {
	  q2 = q1;
	  q1 = q;
	  row2 = row1;
	  row1 = row;
	}
	else if( q > q2 ) {
	  q2 = q;
	  row2 = row;
	}

      } // rows

      double meanrow = sumrow / sumq;
      double rmsrow = 0;
      if( nrow > 1 ) rmsrow = sqrt( sumrow2/sumq - meanrow*meanrow );

      //if( q1 < thrke ) q1 = thrke; // [ke] threshold
      //if( q2 < thrke ) q2 = thrke; // This alone improves resolution !
      if( q1 < 1.5 ) q1 = 1.5; // [ke] threshold like in data
      if( q2 < 1.5 ) q2 = 1.5; // This alone improves resolution !

      double q12 = q1 + q2;

      if( nrow == 2 ) c->row = ( row1*q1 + row2*q2 ) / q12; // overwrite !

      double eta = 0;
      if( q12 > 1 ) eta = ( q1 - q2 ) / q12;
      if( row2 > row1 ) eta = -eta; // for rmsyvseta at 14.5 deg
      //if( row1 > row2 ) eta = -eta; // for rmsyvseta at 19 deg

      if( chip < 0 ) eta = -eta;

      // residuals:
      double dx0 = ( c->col + 0.5 - 0.5*Ncol ) * pitchx - xmid; // [um]
      double dy0 = ( c->row + 0.5 - 0.5*Nrow ) * pitchy - ymid;

      // skew correction
      double skwcorr =  0.0665505 + 77.6218*skw;
      double dx = ( c->col + 0.5 - 0.5*Ncol ) * pitchx - xmid - skwcorr; // [um]
      double dy = ( c->row + 0.5 - 0.5*Nrow ) * pitchy - ymid;

      double dxtr = dx - smr*rx/cos(alfa); // smeared by telescope
      double dytr = dy - smr*ry/cos(beta); // negative sign is important: cancellation with ymod

      double dx0tr = dx0 - smr*rx/cos(alfa); // smeared by telescope
      double dy0tr = dy0 - smr*ry/cos(beta); // negative sign is important: cancellation with ymod

      double dxbg = 180 * ( myRandom->Rndm() - myRandom->Rndm() ); // triangular
      double dybg = 120 * ( myRandom->Rndm() - myRandom->Rndm() ); // [-100,100]

      hdx0.Fill( dx0 );
      hdy0.Fill( dy0 );

      hdx.Fill( dx );
      hdy.Fill( dy );

      if( abs( dy ) < 51 ) hdxc.Fill( dxtr );
      if( abs( dx ) < 76 ) hdyc.Fill( dytr );

      if( qnrm*0.01 > 16 ) {

	if( abs( dy ) < 51 ) hdxcq0.Fill( dxtr );
	if( abs( dx ) < 76 ) hdycq0.Fill( dytr );
	if( qnrm*0.01 < 35 ) {
	  if( abs( dy ) < 51 ) hdxcq1.Fill( dxtr );
	  if( abs( dx ) < 76 ) hdycq1.Fill( dytr );
	}
	if( qnrm*0.01 < 30 ) {
	  if( abs( dy ) < 51 ) hdxcq2.Fill( dxtr );
	  if( abs( dx ) < 76 ) hdycq2.Fill( dytr );
	}
	if( qnrm*0.01 < 25 ) {
	  if( abs( dy ) < 51 ) { 
	    hdxcq3.Fill( dxtr );
	    hdx0cq3.Fill( dx0tr );
	  }
	  if( abs( dx ) < 76 ) {
	    hdycq3.Fill( dytr );
	    hdycq3.Fill( dybg, wbg );
	    if( xmid > -75+40 ) { // no bias dot
	      hdycq3d.Fill( dytr ); // Gaussian smearing
	      hdycq3d.Fill( dybg, wbg );
	    }
	  }
	}

      } // Q0

      if( qnrm*0.01 >  1 && qnrm*0.01 < 25 ) { // 203i
	if( abs( dx ) < 76 ) hdycq5.Fill( dytr );
      }

      // thick 450
      if( qnrm*0.01 > 25 && qnrm*0.01 < 45 ) {
	if( abs( dx ) < 76 ) hdycq450.Fill( dytr );
      }

      // thick 375
      if( qnrm*0.01 > 22 && qnrm*0.01 < 35 ) {
	if( abs( dx ) < 76 ) hdycq375.Fill( dytr );
      }

      // thick 300
      if( qnrm*0.01 > 18 && qnrm*0.01 < 30 ) {
	if( abs( dx ) < 76 ) hdycq300.Fill( dytr );
      }

      // thick 225
      if( qnrm*0.01 > 12 && qnrm*0.01 < 23 ) {
	if( abs( dx ) < 76 ) hdycq225.Fill( dytr );
      }

      // thick 150
      if( qnrm*0.01 >  6 && qnrm*0.01 < 14 ) {
	if( abs( dx ) < 76 ) hdycq150.Fill( dytr );
      }

      // thick  75
      if( qnrm*0.01 >  3 && qnrm*0.01 < 11 ) {
	if( abs( dx ) < 76 ) hdycq075.Fill( dytr );
      }

      // linking cut like in telescope analysis:

      if( abs( dx ) < 150 && abs(dy) < 100 ) {

	hqvis.Fill( qvis*0.01 ); // [ke]
	hqnrm.Fill( qnrm*0.01 ); // normalized Landau
	hqnrmv.Fill( qnrm*0.01 ); // normalized Landau, different binning

	if( xmid > -75+40 ) { // no bias dot
	  hqvis1.Fill( qvis*0.01 ); // [ke]
	  hqvis1v.Fill( qvis*0.01 ); // [ke]
	}

	//if( qnrm*0.01 > 16 && qnrm*0.01 < 30 ) { // lq cut, taken out 5/2014

	qvsxm.Fill( umod+ 75, qvis*0.01 );
	qvsxm.Fill( xmod+225, qvis*0.01 );

	if( xmid > -75+40 ) { // no bias dot
	  qvsym.Fill( ymod+ 50, qvis*0.01 );
	  qvsym.Fill( ymod+150, qvis*0.01 );
	}

	qvsxmym.Fill( umod+ 75, ymod+ 50, qvis*0.01 );
	qvsxmym.Fill( umod+ 75, ymod+150, qvis*0.01 );
	qvsxmym.Fill( xmod+225, ymod+ 50, qvis*0.01 );
	qvsxmym.Fill( xmod+225, ymod+150, qvis*0.01 );

	//} // lq

	hnall.Fill( Nall );
	hnvis.Fill( Nvis );

	npxvsq.Fill( qnrm*0.01, Nvis );
	npxvsqv.Fill( qnrm*0.01, Nvis );

	if( nrow == 2 ) {

	  hq1.Fill(q1);
	  hq2.Fill(q2);
	  rmsyvsq2.Fill( q2, abs(dytr) );
	  rmsyvsq2.Fill( q2, dybg, wbg );
	  q2vsym.Fill( ymod+ 50, q2 );
	  q2vsym.Fill( ymod+150, q2 );

	} // 2-row clus

	rmsxvsq.Fill( qnrm*0.01, abs(dxtr) );
	rmsxvsq.Fill( qnrm*0.01, abs(dxbg), wbg );

	rmsyvsq.Fill( qnrm*0.01, abs(dytr) );
	rmsyvsq.Fill( qnrm*0.01, abs(dybg), wbg );

	rmsyvsqv.Fill( qnrm*0.01, abs(dytr) );
	rmsyvsqv.Fill( qnrm*0.01, abs(dybg), wbg );

	// Landau peak:

	//if( qnrm*0.01 > 18 &&  qnrm*0.01 < 35 ) {
	if( qnrm*0.01 > 16 &&  qnrm*0.01 < 30 ) { // lq cut

	  hncol.Fill( ncol );
	  hnrow.Fill( nrow );
	  
	  // skew
	  skwcol.Fill(skw);
	  skwcolvsxm.Fill( -umod+ 75, skw ); // smeared
	  skwcolvsxm.Fill( xmod+225, skw ); // smeared

	  dxvsskwcol.Fill( skw, dxtr);
	  dx0vsskwcol.Fill( skw, dx0tr);

	  if( xmid > -75+20 && xmid < 75-20 ) { // cut out x-cracks
	    nrowvsy.Fill( ymid+ 50, nrow );
	    nrowvsy.Fill( ymid+150, nrow );
	  }

	  // xmid = -75..75, bias dot not flipped

	  ncolvsxmid.Fill(-xmid+ 75, ncol ); // flipped
	  ncolvsxmid.Fill( xmid+225, ncol );

	  //ncolvsxm.Fill(  75-xmod, ncol ); // smeared
	  ncolvsxm.Fill( umod+ 75, ncol ); // flipped, smeared
	  ncolvsxm.Fill( xmod+225, ncol ); // smeared

	  //nrowvsxm.Fill(  75-xmod, nrow ); // smeared
	  nrowvsxm.Fill( umod+ 75, nrow ); // smeared
	  nrowvsxm.Fill( xmod+225, nrow ); // smeared

	  nrowvsym.Fill( ymod+ 50, nrow );
	  nrowvsym.Fill( ymod+150, nrow );

	  npxvsxmym.Fill( umod+ 75, ymod+ 50, Nvis );
	  npxvsxmym.Fill( umod+ 75, ymod+150, Nvis );
	  npxvsxmym.Fill( xmod+225, ymod+ 50, Nvis );
	  npxvsxmym.Fill( xmod+225, ymod+150, Nvis );

	  if( nrow == 2 ) {
	    heta.Fill( eta );
	    etavsym.Fill(-ymod+ 50, eta );
	    etavsym.Fill(-ymod+150, eta );

	    q0vsym.Fill( ymod+ 50, qrow[rowmin] );
	    q0vsym.Fill( ymod+150, qrow[rowmin] );

	    q9vsym.Fill( ymod+ 50, qrow[rowmax] );
	    q9vsym.Fill( ymod+150, qrow[rowmax] );

	    rmsyvseta.Fill( eta, abs(dytr) );
	    rmsyvseta.Fill( eta, abs(dybg), wbg );
	  }

	  //dyvsxm.Fill(  75-xmod, dytr );
	  dyvsxm.Fill( umod+ 75, dytr );
	  dyvsxm.Fill( xmod+225, dytr );

	  dyvsym.Fill( ymod+ 50, dytr );
	  dyvsym.Fill( ymod+150, dytr );

	  dx0vsxm.Fill( -umod+ 75, dx0tr );
	  dx0vsxm.Fill( xmod+225, dx0tr );

	  dxvsxm.Fill( -umod+ 75, dxtr );
	  dxvsxm.Fill( xmod+225, dxtr );

	  dxvsym.Fill( ymod+ 50, dxtr );
	  dxvsym.Fill( ymod+150, dxtr );

	  //rmsxvsxm.Fill(  75-xmod, abs(dxtr) );
	  rmsxvsxm.Fill( umod+ 75, abs(dxtr) );
	  rmsxvsxm.Fill( umod+ 75, abs(dxbg), wbg );
	  rmsxvsxm.Fill( xmod+225, abs(dxtr) );
	  rmsxvsxm.Fill( xmod+225, abs(dxbg), wbg );

	  //rmsyvsxm.Fill(  75-xmod, abs(dytr) );
	  rmsyvsxm.Fill( umod+ 75, abs(dytr) );
	  rmsyvsxm.Fill( umod+ 75, abs(dybg), wbg );
	  rmsyvsxm.Fill( xmod+225, abs(dytr) );
	  rmsyvsxm.Fill( xmod+225, abs(dybg), wbg );

	  rmsyvsym.Fill(-ymod+ 50, abs(dytr) );
	  rmsyvsym.Fill(-ymod+ 50, abs(dybg), wbg );
	  rmsyvsym.Fill(-ymod+150, abs(dytr) );
	  rmsyvsym.Fill(-ymod+150, abs(dybg), wbg );

	  rmsyvsxmym.Fill( umod+ 75, ymod+ 50, abs(dytr) );
	  rmsyvsxmym.Fill( xmod+225, ymod+ 50, abs(dytr) );
	  rmsyvsxmym.Fill( umod+ 75, ymod+150, abs(dytr) );
	  rmsyvsxmym.Fill( xmod+225, ymod+150, abs(dytr) );

	  if( abs( ymod ) < 5 ) hdymid.Fill( dytr ); // data have resolution dip at mid pix

	} // Landau peak

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ ) {

	  hpxq.Fill( px->anaVcal ); // [ke]
	  if( c->size == 1 ) hpxq1.Fill( px->anaVcal ); // [ke]
	  if( c->size == 2 ) hpxq2.Fill( px->anaVcal ); // [ke]
	  if( c->size >= 3 ) hpxq3.Fill( px->anaVcal ); // [ke]
	  pxqvsq.Fill( qnrm*0.01, px->anaVcal );
	  pxqvsqv.Fill( qnrm*0.01, px->anaVcal );
	  pxqvsxm.Fill( umod+ 75, px->anaVcal );
	  pxqvsxm.Fill( xmod+225, px->anaVcal );
	  pxqvsym.Fill( ymod+ 50, px->anaVcal );
	  pxqvsym.Fill( ymod+150, px->anaVcal );

	} // pixel loop

      } // dx and dy linking cut

    } // clusters

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << endl;
  cout << "eof" << endl;
  cout << "events " << iev << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;
}
