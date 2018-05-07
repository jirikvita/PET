/* File PET.cpp
 *
 * Created       : Wed Oct 31 14:53:32 CDT 2007
 * Author        : kvita
 * Purpose       : 
 * Last modified : 
 * Comments      : 
 */

#include "TSinogram.cpp"
#include "standardsetup.C"

#include <iostream>
#include <string>

#include "TCanvas.h"
//#include "T.h"

using std::endl;
using std::cout;
using std::cerr;

int PET(int nPhiBins = 1,
	int ifilter = 4,
	bool killSinoChannels = false,
	bool smearSino = false,
	double sigma = 0.0005,
	double killValue = 0.0,
	bool noTitle = true)
{
  cout << "=== Running with configuration: ===" << endl
  << " nPhiBins: " << nPhiBins << endl
  << " ifilter: " << ifilter << endl
  << " killSinoChannels: " << killSinoChannels << endl
  << " smearSino: " << smearSino << endl
  << " sigma: " << sigma << endl
  << endl;


  setup();
  std::string option2D = "colz";

  if (noTitle) {
    option2D = "col";
    gStyle -> SetOptTitle(0);
    gStyle->SetPadRightMargin(0.02);//0.02 fot TH1F?!
    gStyle->SetPadLeftMargin(0.05);//0.1
    gStyle->SetPadTopMargin(0.02);//0.1
    gStyle->SetPadBottomMargin(0.05);//0.1
  }

  // must be 2 to the power of 2 4 8 16 or 32!!!
  int nbinsx = 256;
  int nbinsy = 256;
  double x1 = -1., x2 = 1., y1 = -1., y2 = 1.;

  // make the 2D object slice

  // circle+square+ellipse

  TF2 *ellipse = new TF2("Ellipse",
			"[0]*((x-[1])^2/[2]^2+(y-[3])^2/[4]^2 < [5]^2)",
			x1, x2, y1, y2);
  TF2 *square = new TF2("Square",
			"[0]*(x < [1] && x > [2] && y < [3] && y > [4] )",
			x1, x2, y1, y2);
  TF2 *circle = new TF2("Circle",
			"[0]*( ((x-[1])^2 + (y-[2])^2) < [3]^2 )",
			x1, x2, y1, y2);

  //  TF2 *slice1 = new TF2("test_PET_profile", "[0] + Circle[1] + Square(5) + Ellipse(10)", 
  //			x1, x2, y1, y2);

  TF2 *slice1 = new TF2("test_PET_profile",
  			"[0] + "
			"[1]*( ((x-[2])^2 + (y-[3])^2) < [4]^2 ) + "
			"[5]*(x < [6] && x > [7] && y < [8] && y > [9] ) + [10]*((x-[11])^2/[12]^2+(y-[13])^2/[14]^2 < [15]^2) + "
			"[16]*((x-[17])^2/[18]^2+(y-[19])^2/[20]^2 < [21]^2) + "
			"[22]*( ((x-[23])^2 + (y-[24])^2) < [25]^2 ) + "
			"[26]*( ((x-[27])^2 + (y-[28])^2) < [29]^2 )",
  			x1, x2, y1, y2);
  double par[] = {0.01, 
		  1., 0.2, 0.2, 0.3,
		  0.8, -0.2, -0.4, -0.1, -0.3,
		  0.3, 0., 0.9, 0., 0.7, 1.,
		  0.35, 0.23, 0.9, -0.33, 0.90, 0.22,
		  1.1, -0.45, 0.41, 0.03, 
		  1.1, -0.50, 0.482, 0.010 };
  
  for (int ip = 0; ip < slice1 -> GetNpar(); ++ip) { 
    slice1 -> SetParameter(ip, par[ip]);
  }

  int npx = 400;
  int npy = 400;
  slice1 -> SetNpx(npx);
  slice1 -> SetNpy(npy);
   
  TCanvas *can = new TCanvas("PET", "PET", 0, 0, 1000, 800);
  can -> Divide(2, 2);
  TCanvas *can_FFT = new TCanvas("PET_FFT", "PET_FFT", 100, 100, 1000, 800);
  can_FFT -> Divide(2, 2);
  TCanvas *can_FFT_proj = new TCanvas("PET_FFT_proj", "PET_FFT_proj", 200, 200, 1000, 800);
  can_FFT_proj -> Divide(2, 2);



  // the original slice:
  can -> cd(1);
  slice1 -> Draw(option2D.c_str());
  gPad -> Print(Form("can_%s.eps", slice1 -> GetName()));
  gPad -> Print(Form("can_%s.gif", slice1 -> GetName()));


  // make the sinogram
  TSinogram *Sino = new TSinogram(nbinsx, x1, x2, nbinsy, y1, y2, slice1);
  
  const int MaxnPhiBins = 7;
  int PhiBins[MaxnPhiBins] = {64, 
			      2, 4, 8, 
			      16, 32, 128};

  for (int iphi = 0; iphi < nPhiBins; ++iphi) {

    int phibins = PhiBins[iphi];
    int sbins = nbinsx;
    //    double ds = 0.01;
    //    double dz = 0.01;
    Sino -> MakeTheProjection(phibins, sbins);

    if (killSinoChannels) {
      Sino -> KillSinogramChannel(121, killValue);
    }
    if (smearSino) {
      Sino -> SmearSinogram(sigma);
    }

    TH2D *sinogram = Sino -> GetSinogram();
    can -> cd(2);
    sinogram -> SetStats(0);
    sinogram -> Draw(option2D.c_str());
    gPad -> Print(Form("can_%s_%iphi.eps", sinogram -> GetName(), PhiBins[iphi]));
    gPad -> Print(Form("can_%s_%iphi.gif", sinogram -> GetName(), PhiBins[iphi]));

    // make the siple backprojection  
    Sino -> MakeTheBackProjection(Sino -> GetSinogram());
    TH2D *backprojection = Sino -> GetBackProjection();
    can -> cd(3);
    backprojection -> SetStats(0);
    backprojection -> DrawCopy(option2D.c_str());
    gPad -> Print(Form("can_%s_%iphi.eps", backprojection -> GetName(), PhiBins[iphi]));
    gPad -> Print(Form("can_%s_%iphi.gif", backprojection -> GetName(), PhiBins[iphi]));


    TH2D *sinogram_FFT1D_re;
    TH2D *sinogram_FFT1D_im;
    TH2D *sinogram_InvFFT1D_re;
    TH2D *sinogram_InvFFT1D_im;
    TH2D *filtered_backprojection;

    //    if (iphi == 0) {    
    if (1) {

    // 1D FFT of each slice:
    Sino -> MakeSinogramFFT1D();

    sinogram_FFT1D_re = Sino -> GetSinogramFFT1D_re();
    TH1D *sinogram_FFT1D_re_proj;
    if (sinogram_FFT1D_re) {
      can_FFT -> cd(1);
      sinogram_FFT1D_re -> SetStats(0);
      sinogram_FFT1D_re -> Draw(option2D.c_str());
      gPad -> Print(Form("can_%s_%iphi.eps", sinogram_FFT1D_re -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi.gif", sinogram_FFT1D_re -> GetName(), PhiBins[iphi]));

      int N = sinogram_FFT1D_re -> GetNbinsX();
      sinogram_FFT1D_re -> GetXaxis() -> SetRange(0, N/8);
      gPad -> Print(Form("can_%s_%iphi_low.eps", sinogram_FFT1D_re -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi_low.gif", sinogram_FFT1D_re -> GetName(), PhiBins[iphi]));
      sinogram_FFT1D_re -> GetXaxis() -> SetRange(7*N/8, N);
      gPad -> Print(Form("can_%s_%iphi_high.eps", sinogram_FFT1D_re -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi_high.gif", sinogram_FFT1D_re -> GetName(), PhiBins[iphi]));
      sinogram_FFT1D_re -> GetXaxis() -> SetRange(3*N/8, 5*N/8);
      gPad -> Print(Form("can_%s_%iphi_middle.eps", sinogram_FFT1D_re -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi_middle.gif", sinogram_FFT1D_re -> GetName(), PhiBins[iphi]));

      can_FFT_proj -> cd(1);
      sinogram_FFT1D_re_proj =  (TH1D*) sinogram_FFT1D_re -> ProjectionX("proj_re");
      sinogram_FFT1D_re_proj -> SetStats(0);
      sinogram_FFT1D_re_proj -> Draw("hist");
    } else
      cerr << "ERROR getting sinogram_FFT1D_re!" << endl;
    
    sinogram_FFT1D_im = Sino -> GetSinogramFFT1D_im();
    TH1D *sinogram_FFT1D_im_proj;
    if (sinogram_FFT1D_im) {
      can_FFT -> cd(2);
      sinogram_FFT1D_im -> SetStats(0);
      sinogram_FFT1D_im -> Draw(option2D.c_str());
      gPad -> Print(Form("can_%s_%iphi.eps", sinogram_FFT1D_im -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi.gif", sinogram_FFT1D_im -> GetName(), PhiBins[iphi]));

      int N = sinogram_FFT1D_im -> GetNbinsX();
      sinogram_FFT1D_im -> GetXaxis() -> SetRange(0, N/8);
      gPad -> Print(Form("can_%s_%iphi_low.eps", sinogram_FFT1D_im -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi_low.gif", sinogram_FFT1D_im -> GetName(), PhiBins[iphi]));
      sinogram_FFT1D_im -> GetXaxis() -> SetRange(7*N/8, N);
      gPad -> Print(Form("can_%s_%iphi_high.eps", sinogram_FFT1D_im -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi_high.gif", sinogram_FFT1D_im -> GetName(), PhiBins[iphi]));
      sinogram_FFT1D_im -> GetXaxis() -> SetRange(3*N/8, 5*N/8);
      gPad -> Print(Form("can_%s_%iphi_middle.eps", sinogram_FFT1D_im -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi_middle.gif", sinogram_FFT1D_im -> GetName(), PhiBins[iphi]));

      can_FFT_proj -> cd(2);
      sinogram_FFT1D_im_proj =  (TH1D*) sinogram_FFT1D_im -> ProjectionX("proj_im");
      sinogram_FFT1D_im_proj -> SetStats(0);
      sinogram_FFT1D_im_proj -> Draw("hist");

    } else
      cerr << "ERROR getting sinogram_FFT1D_im!" << endl;


    // Apply filter !!!
    Sino -> ApplyFilter(ifilter);

    // Inverse 1D FFT of each slice:
    Sino -> MakeSinogramInvFFT1D();

    can_FFT -> cd(3);
    sinogram_InvFFT1D_re = Sino -> GetSinogramInvFFT1D_re();
    if (sinogram_InvFFT1D_re) {
      sinogram_InvFFT1D_re -> SetStats(0);
      sinogram_InvFFT1D_re -> Draw(option2D.c_str());
      gPad -> Print(Form("can_%s_%iphi.eps", sinogram_InvFFT1D_re -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi.gif", sinogram_InvFFT1D_re -> GetName(), PhiBins[iphi]));
    } else
      cerr << "ERROR getting sinogram_InvFFT1D_re!" << endl;
    
    can_FFT -> cd(4);
    sinogram_InvFFT1D_im = Sino -> GetSinogramInvFFT1D_im();
    if (sinogram_InvFFT1D_im) {
      sinogram_InvFFT1D_im -> SetStats(0);
      sinogram_InvFFT1D_im -> Draw(option2D.c_str());
      gPad -> Print(Form("can_%s_%iphi.eps", sinogram_InvFFT1D_im -> GetName(), PhiBins[iphi]));
      gPad -> Print(Form("can_%s_%iphi.gif", sinogram_InvFFT1D_im -> GetName(), PhiBins[iphi]));
    } else
      cerr << "ERROR getting sinogram_InvFFT1D_im!" << endl;


    // make backprojection with the inversed image:
    Sino -> MakeTheBackProjection(Sino -> GetSinogramInvFFT1D_re());
    filtered_backprojection = Sino -> GetBackProjection();
    can -> cd(4);
    filtered_backprojection -> SetStats(0);
    filtered_backprojection -> DrawCopy(option2D.c_str());
    gPad -> Print(Form("can_%s_%iphi.eps", filtered_backprojection -> GetName(), PhiBins[iphi]));
    gPad -> Print(Form("can_%s_%iphi.gif", filtered_backprojection -> GetName(), PhiBins[iphi]));



    // smear, blur...

    // deconvolute
    // (using filters)

    // compare to the original 2D slice


    } // iphi == 0


    can -> Print(Form("%s_%iiphi.eps", can -> GetName(), PhiBins[iphi]));
    can -> Print(Form("%s_%iiphi.gif", can -> GetName(), PhiBins[iphi]));

    can_FFT -> Print(Form("%s_%iiphi.eps", can_FFT -> GetName(), PhiBins[iphi]));
    can_FFT -> Print(Form("%s_%iiphi.gif", can_FFT -> GetName(), PhiBins[iphi]));

    can_FFT_proj -> Print(Form("%s_%iiphi.eps", can_FFT_proj -> GetName(), PhiBins[iphi]));
    can_FFT_proj -> Print(Form("%s_%iiphi.gif", can_FFT_proj -> GetName(), PhiBins[iphi]));

  } // phi bins
  int MaxN =  Sino -> GetSinogram() -> GetNbinsX();
  TH1D *filter_h = Sino -> MakeFilterHisto(ifilter, MaxN, Form("Filter%i", ifilter), Form("Filter%i", ifilter));
  can_FFT_proj -> cd(3);
  filter_h -> SetStats(0);
  filter_h -> Draw("hist");
  gPad -> Print(Form("filter_%s.eps", filter_h -> GetName()));
  gPad -> Print(Form("filter_%s.gif", filter_h -> GetName()));


  return 0;

}
