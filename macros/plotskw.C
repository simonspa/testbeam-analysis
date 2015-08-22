void plotskw() {

  TCanvas *c1 = new TCanvas("c1","par0",600,600);
  TProfile *ppar0 = new TProfile("par0"," ",90,0,90,-10,10,"");

  TCanvas *c2 = new TCanvas("c2","par1",600,600);
  TProfile *ppar1 = new TProfile("par1"," ",90,0,90,-2000,1000,"");


  Double_t tilt, par0, par1;

  ifstream in;
  string line;
  in.open("simulation/skwcorr_504_294um_gap30.dat");

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
    ppar0->Fill(tilt,par0,1);
    ppar1->Fill(tilt,par1,1);
  }


  c1->cd();
  ppar0->Draw();

  c2->cd();
  ppar1->Draw();

  in.close();
}
