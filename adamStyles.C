/////////
// the Style Section
/////////
void setHLTStyle() {
  setTDRStyle();
  //  TStyle *hltStyle = new TStyle(*tdrStyle);
  TStyle *hltStyle  = new TStyle("hltStyle","My HLT Styles");
  gStyle->Copy(*hltStyle);

  hltStyle->SetCanvasColor(-1);

  hltStyle->SetCanvasDefH(600);
  hltStyle->SetCanvasDefW(600);
  hltStyle->SetPadColor(-1);  
  hltStyle->SetPadGridX(false);
  hltStyle->SetPadGridY(false);
  hltStyle->SetGridWidth(0.25);
  
  hltStyle->SetFrameFillColor(-1); // Transparent

  hltStyle->SetHistFillColor(-1); // Transparent
  hltStyle->SetHistFillStyle(0);  // None
  hltStyle->SetHistLineWidth(3);

  hltStyle->SetPadTopMargin(0.08);
  hltStyle->SetPadBottomMargin(0.12);
  hltStyle->SetPadLeftMargin(0.15);
  hltStyle->SetPadRightMargin(0.04);

  hltStyle->SetTitleSize(0.05);
  hltStyle->SetTitleFillColor(-1); // Transparent
  hltStyle->SetTitleH(0.05); // Set the height of the title box
  hltStyle->SetTitleW(0.); // Set the width of the title box

  hltStyle->SetTitleSize(0.04, "XYZ");     
  hltStyle->SetTitleOffset(1.2, "X"); // Another way to set the Offset
  hltStyle->SetTitleOffset(1.8, "Y"); // Another way to set the Offset

  hltStyle->SetLabelSize(0.035, "XYZ");

  hltStyle->SetPalette(1,0);
  hltStyle->SetFillColor(0);    // White
  hltStyle->SetFillStyle(4000); // Transparent

  hltStyle->SetStatStyle(0);
  hltStyle->SetTitleStyle(0);
  hltStyle->SetCanvasBorderSize(0);
  hltStyle->SetFrameBorderSize(0);
  hltStyle->SetLegendBorderSize(0);
  hltStyle->SetStatBorderSize(0);
  hltStyle->SetTitleBorderSize(0);

  hltStyle->cd();
}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(500); //Height of canvas for form=1 //600
  tdrStyle->SetCanvasDefW(700); //Width of canvas  for form=1 //600
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(true);
  tdrStyle->SetPadGridY(true); //grid y axis
//  tdrStyle->SetGridColor(0);
//  tdrStyle->SetGridStyle(3);
//  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  //aaa tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  //aaa  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  //alexey
  //tdrStyle->SetOptFit(1);
  tdrStyle->SetOptFit(1111);
  tdrStyle->SetFitFormat("5.4g");
  //aaa tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(1); // To display the mean and RMS:   SetOptStat("mr");
  //alexey 1 line
  //tdrStyle->SetOptStat(0110);
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.045);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.15); //aaa 0.13
  tdrStyle->SetPadRightMargin(0.03); //aaa 0.05

// For the Global title:

   //alexey 1 line
  tdrStyle->SetOptTitle(kTRUE);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.17); //aaa 1.05
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(true);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(false); //right ruler yes/no

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}

void setRBellanStyle() {
  TStyle *theStyle = new TStyle("rbStyle", "Style for Bellan Thesis");
  theStyle->SetOptStat(1);
  //theStyle->SetOptStat(0110);
  theStyle->SetPadBorderMode(0);
  theStyle->SetCanvasBorderMode(0);
  theStyle->SetPadColor(0);
  theStyle->SetCanvasColor(0);
  theStyle->SetMarkerStyle(8);
  theStyle->SetMarkerSize(0.7);
  theStyle->SetPalette(1);
  
  theStyle->SetStatH(0.3);
  //   theStyle->SetTextFont(132);
  //   theStyle->SetTitleFont(132);
  theStyle->SetTitleBorderSize(1);
  //    theStyle->SetPalette(1);
  theStyle->SetOptStat(1);
  //theStyle->SetOptStat(0110);
  theStyle->SetFitFormat("4.4g");
  theStyle->SetStatY(0.99);
  theStyle->SetStatX(0.99);
  theStyle->SetTitleYOffset(1.6);
  theStyle->SetLabelSize(0.035, "XYZ");
  theStyle->SetPadGridX(true);
  theStyle->SetPadGridY(true);
  theStyle->SetFrameBorderMode(0);
  theStyle->SetTitleFillColor(0);
  theStyle->SetLegendBorderSize();
  
  // theStyle->SetCanvasDefH(600);
  // theStyle->SetCanvasDefW(400);
  
  //theStyle->SetOptLogy(); //aaa
  // theStyle->SetOptLogx();
  theStyle->cd();
}
