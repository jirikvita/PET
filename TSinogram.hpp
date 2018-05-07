/* File          : TSinogram.hpp
 *
 * Created       : Wed Oct 31 14:53:32 CDT 2007
 * Author        : kvita
 * Purpose       : 
 * Last modified : 
 * Comments      : 
 */

#ifndef TSINOGRAM_HPP
#define TSINOGRAM_HPP

#include "TH2D.h"
#include "TH1D.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom3.h"
// #include "T.h"


#include "macroDensity.C"

#include <iostream>

#include "TCanvas.h"
//#include "T.h"

using std::endl;
using std::cout;
using std::cerr;


const int maxNbins = 1024;

class TSinogram {

public:
  
  TSinogram(int binsx, double x1, double x2, int binsy, double y1, double y2);
  TSinogram(int binsx, double x1, double x2, int binsy, double y1, double y2, TF2 *slice);
  void Init(int binsx, double x1, double x2, int binsy, double y1, double y2, TF2 *slice);
  ~TSinogram();

  // methods:
  void MakeTheProjection(int phibins, int sbins, double dzSF = 1., int nsubsbins = 1);
  void MakeTheBackProjection(TH2D *sinogram);
  void MakeSinogramFFT1D();
  void MakeSinogramInvFFT1D();

  double EvalFilter(int ifilter, int channel, int maxN);
  TH1D* MakeFilterHisto(int ifilter, int maxN, std::string name, std::string title);
  void ApplyFilter(int ifilter);
  void KillSinogramChannel(int channel, double value);
  void SmearSinogram(double sigma);
  void Reconstruct(); // perform the deconvolution
  void DrawSinogram();
  void DrawReconstructed();

  bool CheckPower(int FTN);
  int GetReversedBitOrder(int m, int power);
  int CheckPowerBinning(int base, int n);
  
  // setters:

  void SetOriginalSlice(TF2 *slice);

  // getters:

  TH2D* GetSinogram();
  TH2D* GetBackProjection();
  TH2D* GetSinogramFFT1D_re();
  TH2D* GetSinogramFFT1D_im();

  TH2D* GetSinogramInvFFT1D_re();
  TH2D* GetSinogramInvFFT1D_im();

  TF2* GetOriginalSlice();
  TH2D* GetBinnedSlice();
  TH2D* GetReconstructedImage();

  // try to unfold the smearing later?

  // Simulate smearing due to:
  // - attenuation bias
  // - finite positron path
  // - collinearity bias
  // - noise
  // - coincidence fakes

private:

  int _binsx;
  int _binsy;
  double _x1, _x2, _y1, _y2;
  double _sigma; // for smearing

  TRandom3 *_rand;

  // data members
  TF2 *_slice_f2;
  TH2D *_binned_slice_h2; // a discrete version
  TH2D *_reconstructed;
  TH2D *_sinogram;
  TH2D *_backprojection;

  // individual slices:
  TH1D *_sinogram_FFT1D_re_h1[maxNbins];
  TH1D *_sinogram_FFT1D_im_h1[maxNbins];

  // sinogram 1D FFT slices in 2D histo:
  TH2D *_sinogram_FFT1D_re;
  TH2D *_sinogram_FFT1D_im;

  double *_sinogram_re[maxNbins];
  double *_sinogram_im[maxNbins];

  // inverse FFT:
  TH2D *_sinogram_InvFFT1D_re;
  TH2D *_sinogram_InvFFT1D_im;

  double *_sinogram_Inv_re[maxNbins];
  double *_sinogram_Inv_im[maxNbins];
  

};


#endif
