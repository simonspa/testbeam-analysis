#include <TPaveText.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TLegend.h>
#include <TMath.h>

void setHHStyle(TStyle& HHStyle)
{
    const int fontstyle=42;
    HHStyle.SetPalette(1);
        
    // ==============
    //  Canvas
    // ==============
            
    HHStyle.SetCanvasBorderMode(0);
    HHStyle.SetCanvasColor(kWhite);
    HHStyle.SetCanvasDefH(600); //Height of canvas
    HHStyle.SetCanvasDefW(600); //Width of canvas
    HHStyle.SetCanvasDefX(0);   //Position on screen
    HHStyle.SetCanvasDefY(0);
            
    // ==============
    //  Pad
    // ==============
            
    HHStyle.SetPadBorderMode(0);
    // HHStyle.SetPadBorderSize(Width_t size = 1);
    HHStyle.SetPadColor(kWhite);
    HHStyle.SetPadGridX(false);
    HHStyle.SetPadGridY(false);
    HHStyle.SetGridColor(kGray);
    HHStyle.SetGridStyle(3);
    HHStyle.SetGridWidth(1);
            
    // ==============
    //  Frame
    // ==============
            
    HHStyle.SetFrameBorderMode(0);
    HHStyle.SetFrameBorderSize(1);
    HHStyle.SetFrameFillColor(0);
    HHStyle.SetFrameFillStyle(0);
    HHStyle.SetFrameLineColor(1);
    HHStyle.SetFrameLineStyle(1);
    HHStyle.SetFrameLineWidth(1);
            
    // ==============
    //  Histo
    // ==============

    HHStyle.SetErrorX(0.0);
    HHStyle.SetEndErrorSize(8);
            
    // HHStyle.SetHistFillColor(1);
    // HHStyle.SetHistFillStyle(0);
    // HHStyle.SetHistLineColor(1);
    HHStyle.SetHistLineStyle(0);
    HHStyle.SetHistLineWidth(1);
    // HHStyle.SetLegoInnerR(Float_t rad = 0.5);
    // HHStyle.SetNumberContours(Int_t number = 20);

    // HHStyle.SetErrorMarker(20);
            
    HHStyle.SetMarkerStyle(20);
            
    // ==============
    //  Fit/function
    // ==============
            
    HHStyle.SetOptFit(0);
    HHStyle.SetFitFormat("5.4g");
    HHStyle.SetFuncColor(2);
    HHStyle.SetFuncStyle(1);
    HHStyle.SetFuncWidth(1);
            
    // ==============
    //  Date
    // ============== 
            
    HHStyle.SetOptDate(0);
    // HHStyle.SetDateX(Float_t x = 0.01);
    // HHStyle.SetDateY(Float_t y = 0.01);
            
    // =====================
    //  Statistics Box
    // =====================
            
    HHStyle.SetOptFile(0);
    HHStyle.SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    HHStyle.SetStatColor(kWhite);
    HHStyle.SetStatFont(fontstyle);
    HHStyle.SetStatFontSize(0.025);
    HHStyle.SetStatTextColor(1);
    HHStyle.SetStatFormat("6.4g");
    HHStyle.SetStatBorderSize(1);
    HHStyle.SetStatH(0.1);
    HHStyle.SetStatW(0.15);
    // HHStyle.SetStatStyle(Style_t style = 1001);
    // HHStyle.SetStatX(Float_t x = 0);
    // HHStyle.SetStatY(Float_t y = 0);
            
    // ==============
    //  Margins
    // ==============

    HHStyle.SetPadTopMargin(0.1);
    HHStyle.SetPadBottomMargin(0.15);
    HHStyle.SetPadLeftMargin(0.20);
    HHStyle.SetPadRightMargin(0.05);
            
    // ==============
    //  Global Title
    // ==============
            
    HHStyle.SetOptTitle(0);
    HHStyle.SetTitleFont(fontstyle);
    HHStyle.SetTitleColor(1);
    HHStyle.SetTitleTextColor(1);
    HHStyle.SetTitleFillColor(10);
    HHStyle.SetTitleFontSize(0.05);
    // HHStyle.SetTitleH(0); // Set the height of the title box
    // HHStyle.SetTitleW(0); // Set the width of the title box
    // HHStyle.SetTitleX(0); // Set the position of the title box
    // HHStyle.SetTitleY(0.985); // Set the position of the title box
    // HHStyle.SetTitleStyle(Style_t style = 1001);
    // HHStyle.SetTitleBorderSize(2);
            
    // ==============
    //  Axis titles
    // ==============
            
    HHStyle.SetTitleColor(1, "XYZ");
    HHStyle.SetTitleFont(fontstyle, "XYZ");
    HHStyle.SetTitleSize(0.04, "XYZ");
    // HHStyle.SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // HHStyle.SetTitleYSize(Float_t size = 0.02);
    HHStyle.SetTitleXOffset(1.0);
    HHStyle.SetTitleYOffset(1.7);
    // HHStyle.SetTitleOffset(1.1, "Y"); // Another way to set the Offset
            
    // ==============
    //  Axis Label
    // ==============
            
    //HHStyle.SetLabelColor(1, "XYZ");
    HHStyle.SetLabelFont(fontstyle, "XYZ");
    HHStyle.SetLabelOffset(0.007, "XYZ");
    HHStyle.SetLabelSize(0.04, "XYZ");
            
    // ==============
    //  Axis
    // ==============
            
    HHStyle.SetAxisColor(1, "XYZ");
    HHStyle.SetStripDecimals(kTRUE);
    HHStyle.SetTickLength(0.03, "XYZ");
    HHStyle.SetNdivisions(510, "XYZ");
    HHStyle.SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    HHStyle.SetPadTickY(1);
            
    // Change for log plots:
    HHStyle.SetOptLogx(0);
    HHStyle.SetOptLogy(0);
    HHStyle.SetOptLogz(0);
            
    // ==============
    //  Text
    // ==============
            
    HHStyle.SetTextAlign(11);
    HHStyle.SetTextAngle(0);
    HHStyle.SetTextColor(1);
    HHStyle.SetTextFont(fontstyle);
    HHStyle.SetTextSize(0.04);
            
    // =====================
    //  Postscript options:
    // =====================
            
    HHStyle.SetPaperSize(20.,20.);
    // HHStyle.SetLineScalePS(Float_t scale = 3);
    // HHStyle.SetLineStyleString(Int_t i, const char* text);
    // HHStyle.SetHeaderPS(const char* header);
    // HHStyle.SetTitlePS(const char* pstitle);
            
    // HHStyle.SetBarOffset(Float_t baroff = 0.5);
    // HHStyle.SetBarWidth(Float_t barwidth = 0.5);
    // HHStyle.SetPaintTextFormat(const char* format = "g");
    // HHStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // HHStyle.SetTimeOffset(Double_t toffset);
    // HHStyle.SetHistMinimumZero(kTRUE);
}

// Draw label for Decay Channel in upper left corner of plot
void DrawPrelimLabel(int cmsprelim, double textSize) {

    TPaveText *cms = new TPaveText();
    cms->AddText("CMS");

    cms->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    cms->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    cms->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    cms->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    cms->SetFillStyle(0);
    cms->SetBorderSize(0);
    if (textSize!=0) cms->SetTextSize(textSize*1.1);
    cms->SetTextAlign(12);
    cms->SetTextFont(61);
    cms->Draw("same");

    if(cmsprelim > 0) {
      TPaveText *extra = new TPaveText();
      if(cmsprelim == 2) { extra->AddText("Private Work"); }
      else { extra->AddText("Preliminary"); }

      extra->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
      extra->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.10 );
      extra->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
      extra->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );

      extra->SetFillStyle(0);
      extra->SetBorderSize(0);
      if (textSize!=0) extra->SetTextSize(textSize);
      extra->SetTextAlign(12);
      extra->SetTextFont(52);
      extra->Draw("same");
    }
}

// Draw statistics info above plot
void DrawCMSLabels(double lumi, double energy, double textSize) {

  //const char *text = "%2.1f #times 10^{6} clusters (fiducial) (%1.1f GeV)";
  const char *text;
  if(lumi > 10000000) {
    lumi /= 10000000;
    text = "%2.1f #times 10^{7} tracks (%1.1f GeV)";
  }
  else {
    lumi /= 1000000;
    text = "%2.1f #times 10^{6} tracks (%1.1f GeV)";
  }
    
  TPaveText *label = new TPaveText();
  label->SetX1NDC(gStyle->GetPadLeftMargin());
  label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
  label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
  label->SetY2NDC(1.0);
  label->SetTextFont(42);
  label->AddText(Form(text, lumi, energy));
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  if (textSize!=0) label->SetTextSize(textSize);
  label->SetTextAlign(32);
  label->Draw("same");
}

void DrawFreeCMSLabels(char* text, double textSize) {

  //const char *text = "%2.1f #times 10^{6} clusters (fiducial) (%1.1f GeV)";
  //char *text = "%2.1f #times 10^{6} clusters (fiducial)";
    
  TPaveText *label = new TPaveText();
  label->SetX1NDC(gStyle->GetPadLeftMargin());
  label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
  label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
  label->SetY2NDC(1.0);
  label->SetTextFont(42);
  label->AddText(text);
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  if (textSize!=0) label->SetTextSize(textSize);
  label->SetTextAlign(32);
  label->Draw("same");
}

void setLegendStyle(TLegend *leg) {

  Double_t x1 = 0.5216107;
  Double_t y1 = 0.7839721;
  double height = 0.075, width = 0.275;

  leg->SetX1NDC(x1);
  leg->SetY1NDC(y1);
  leg->SetX2NDC(x1 + width);
  leg->SetY2NDC(y1 + height);

  leg->SetTextFont(42);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
}

void setStyle(TProfile *hist, TString name)
{
  hist->SetLineWidth(1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.08);
  hist->GetYaxis()->SetTitleOffset(1.7);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelOffset(0.007);

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->SetLineWidth(3);
  
  if(name == "data") {
    hist->SetMarkerColor(kBlack);
    hist->SetFillStyle(3005);
    hist->SetFillColor(kBlack);
  }
  else if(name == "dot1") {
    hist->SetLineWidth(3);
    hist->SetMarkerSize(0);
    hist->SetLineColor(kCyan+1);
    hist->SetFillStyle(3004);
    hist->SetFillColor(kRed+1);
  }
  else {
    hist->SetLineWidth(3);
    hist->SetMarkerSize(0);
    hist->SetLineColor(kRed+1);
    hist->SetFillStyle(3004);
    hist->SetFillColor(kRed+1);
  }
}

void setStyle(TGraph *hist, TString name)
{
  hist->SetLineWidth(1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.08);
  hist->GetYaxis()->SetTitleOffset(1.7);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelOffset(0.007);

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->SetLineWidth(3);

  if(name == "data") {
    hist->SetMarkerColor(kBlack);
    hist->SetFillStyle(3005);
    hist->SetFillColor(kBlack);
  }
  else if(name == "dot1") {
    hist->SetLineWidth(3);
    hist->SetMarkerSize(0);
    hist->SetLineColor(kCyan+1);
    hist->SetFillStyle(3004);
    hist->SetFillColor(kRed+1);
  }
  else {
    hist->SetLineWidth(3);
    hist->SetMarkerSize(0);
    hist->SetLineColor(kRed+1);
    hist->SetFillStyle(3004);
    hist->SetFillColor(kRed+1);
  }
}

void setStyleAndFillLegend(TProfile* hist, TString name, TLegend *leg) {

  setStyle(hist,name);

  if(name == "data") { if(leg) leg->AddEntry(hist, "Data",  "p"); }
  else if(name == "dot1") { if(leg) leg->AddEntry(hist, "Simulation (pixelav) - dot1 design",  "l"); }
  else if(name == "gap30") { if(leg) leg->AddEntry(hist, "Simulation (pixelav) - gap30 design",  "l"); }
  else { if(leg) leg->AddEntry(hist, "Simulation (pixelav)",  "l"); }
}

void setStyleAndFillLegend(TGraph* hist, TString name, TLegend *leg) {

  setStyle(hist,name);

  if(name == "data") { if(leg) leg->AddEntry(hist, "Data",  "p"); }
  else if(name == "dot1") { if(leg) leg->AddEntry(hist, "Simulation (pixelav) - dot1 design",  "l"); }
  else if(name == "gap30") { if(leg) leg->AddEntry(hist, "Simulation (pixelav) - gap30 design",  "l"); }
  else { if(leg) leg->AddEntry(hist, "Simulation (pixelav)",  "l"); }
}
