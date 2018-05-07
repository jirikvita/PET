


#include <iostream>

#include "macroDensity.C"
#include "standardsetup.C"

using std::endl;
using std::cout;
using std::cerr;

void test()

{


  setup();

  double x1 = -2*TMath::Pi();
  double x2 = 2*TMath::Pi();

  TF1 *f = new TF1("testFun", "[0] + [1]*sin([2]*x + [3]) + [4]*sin([5]*x + [6]) + [7]*sin([8]*x + [9])", x1, x2);
  f -> SetParameters(0.,
		     1., 1., 0., 
		     1., 2, 0.,
		     1.5, 3., 0.);
 
  int FTN = 8;
  int nbins = 1;
  for (int j = 0; j < FTN; ++j)
    nbins *= 2;

  TH1D *h = new TH1D("testHisto", "testHiosto", nbins, x1, x2);	    
  //  h -> FillRandom(f -> GetName());
   for (int i = 1; i <= nbins; ++i) {
     h -> SetBinContent(i, f -> Eval(h -> GetBinCenter(i)));
   }

  cout << "Preparing vectors..." << endl;
  double *re = new double[nbins];
  // h -> GetYaxis() -> GetXbins() -> fArray;
  //  int n = h -> GetXaxis() -> GetNbins();
  double *im = new double[nbins];
  for (int i = 0; i < nbins; ++i) {
    re[i] = h -> GetBinContent(i+1);
    im[i] = 0.;
  }

  // see macroDensity.C:
  //  FFT(Int_t n, Double_t *f_re, Double_t *f_im)
  cout << "Running FFT..." << endl;
  FFT(FTN, re, im);
  cout << "Done!" << endl;


  TH1D *fft_re = new TH1D("fft_re", "fft_re", nbins, 0, nbins);
  TH1D *fft_im = new TH1D("fft_im", "fft_im", nbins, 0, nbins);
  for (int i = 0; i < nbins; ++i) {
    fft_re -> SetBinContent(i, re[i]);
    fft_im -> SetBinContent(i, im[i]);
    cout << "re: " << re[i] << " im: " << im[i] << endl;
  }

  TCanvas *can = new TCanvas("FFT", "FFT", 0, 0, 1000, 800);
  can -> Divide(2, 2);
  
  can -> cd(1);
  f -> Draw();

  can -> cd(2);
  h -> Draw("hist");

  can -> cd(3);
  fft_re -> Draw("hist");
  can -> cd(4);
  fft_im -> Draw("hist");

  
  
  

}
