// jiri kvita 2004-2007
// kvita@fnal.gov

#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"

/*
int kOrange;
int kLightOrange;
int kDarkOrange;
int kGreenBlue;
int kLightBlue;
int  kGreenYellow;
int kPurple;
int kDarkGreen;
int kDarkBlue;
int kDarkRed;
int kKhaki;
int kDarkPurple;
*/

void setup()
{

  /*
  TColor colorClass;
  kOrange       = colorClass.GetColor("#ee6f00");
  kLightOrange  = colorClass.GetColor("#ffaa00");
  kDarkOrange   = colorClass.GetColor("#ff5a00");
  kGreenBlue    = colorClass.GetColor("#00c081");
  kLightBlue    = colorClass.GetColor("#00a8ff");
  kGreenYellow  = colorClass.GetColor("#a8ff00");
  kPurple       = colorClass.GetColor("#9c00ff");
  kDarkGreen    = colorClass.GetColor("#008f00");
  kDarkBlue     = colorClass.GetColor("#000085");
  kDarkRed      = colorClass.GetColor("#9f0000");
  kKhaki        = colorClass.GetColor("#848400");
  kDarkPurple   = colorClass.GetColor("#450084");
  */


 TStyle *plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");
  
 //some methods doesn't work in plain, but work called for gStyle...happy ROOTing!!!

  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetTitleColor(0);
  plain->SetStatColor(0);
  plain->SetOptFit(1111);

 
  //  plain->SetHistLineWidth(2);
  //  plain->SetHistLineStyle(0);
  //  plain->SetStatX(0.1);
  
  //plain->SetTitleXSize(0.06);
  //plain->SetTitleYSize(0.06);
  //plain->SetTitleFontSize(0.06);
  //plain->SetTitleOffset(0.5);

  //  plain->SetStatW(0.4);//0.19
  //plain->SetStatH(0.5);//0.1

   //plain->SetPadGridX(true);
  //plain->SetPadGridY(true);

  gROOT->SetStyle("Plain");


  //gStyle->SetOptStat(111111);//all but bins integral
  //gStyle->SetOptStat(1111);//entries,mean,rms, name
  gStyle->SetOptStat(1110);//entries,mean,rms, NO name!
  
  //gStyle->SetOptFit(1100);

  //gStyle ->SetStatFontSize(0.03);
  gStyle->SetStatW(0.19);//0.19
  gStyle->SetStatH(0.2);//0.25
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.9);//0.99
  gStyle->SetStatFontSize(0.08);

  gStyle->SetPadRightMargin(0.12);//0.02 fot TH1F?!
  gStyle->SetPadLeftMargin(0.15);//0.1
  gStyle->SetPadBottomMargin(0.16);//0.1

  // !!! aaa
  //  gStyle->SetPadGridX(true);
  // gStyle->SetPadGridY(true);

  //histogram axs titles
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleYSize(0.06);
 
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleOffset(0.8,"X");
  gStyle->SetTitleOffset(0.8,"Y");
  
  //title of the histogram
  gStyle->SetTitleH(0.09);
  gStyle->SetTitleW(0.73);//0.6
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleY(0.99);
  gStyle->SetPalette(1);
  gStyle->SetFillColor(10);

  gROOT->ForceStyle();
}
