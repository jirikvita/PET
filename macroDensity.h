/******************************************************************************
* macro to test median polish                                                 *
*                                                                             *
* Author: Dr. Christian Stratowa, Vienna, Austria.                            *
* Created: 02 Oct 2004                             Last modified: 03 Oct 2004 *
*
* take from file: weightedkerneldensity.c
* Copyright (C) 2003 Ben Bolstad
*
* aim : compute weighted kernel density estimates
*
* the aim here is to implement kernel density estimators, with the option to 
* weight each observation. For speed we will use the FFT to convolve a weighted
* histogram with a kernel.
* 
******************************************************************************/

#ifndef MACRODENSITY_H
#define MACRODENSITY_H

/* 

jiri's notes:

Example running:
   .L macroDensity.C
   TestDensity("epanechnikov") 

taken from
http://root.cern.ch/phpBB2/viewtopic.php?t=2462&sid=5b3abf36e0fadbcbb09cea99e6fc9170

created using
http://stat.ethz.ch/R-manual/R-patched/library/stats/html/fft.html
http://stat.ethz.ch/R-manual/R-patched/library/stats/html/density.html

*/


#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"


#include <Riostream.h>


//______________________________________________________________________________
//Double_t TStat::Mean(Int_t n, const Double_t *arr)
Double_t Mean(Int_t n, const Double_t *arr);
//______________________________________________________________________________
//Double_t TStat::Var(Int_t n, const Double_t *arr, const Double_t mean)
Double_t Var(Int_t n, const Double_t *arr, const Double_t mean);
//______________________________________________________________________________
void MassDist(Int_t nx, Double_t *x, Double_t *w, Double_t xmin, Double_t xmax,
              Int_t ny, Double_t *y);

//______________________________________________________________________________
void TwiddleFactor4FFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im);
//______________________________________________________________________________
void TwiddleFactor4IFFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im);
//______________________________________________________________________________
void FFT(Int_t n, Double_t *f_re, Double_t *f_im);
//______________________________________________________________________________
void IFFT(Int_t n, Double_t *f_re, Double_t *f_im);

//______________________________________________________________________________
Int_t FFTDensityConvolve(Int_t n, Double_t *x_re, Double_t *y_re);

//______________________________________________________________________________
void Kernelize(Int_t n, Double_t *x, Double_t bw, const char *kernel);
//______________________________________________________________________________
Double_t Bandwidth(Int_t n, Double_t *x, Double_t iqr);
//______________________________________________________________________________
void LinearInterpolate(Double_t *xin, Double_t *yin,
                       Int_t nout, Double_t *xout, Double_t *yout);
//______________________________________________________________________________
Int_t Density(Int_t n, Double_t *x, Double_t *w, Int_t nout,  Double_t *xout,
              Double_t *yout, const char *kernel = "epanechnikov");
//______________________________________________________________________________
void TestDensity(const char *kernel = "epanechnikov");

//______________________________________________________________________________
void TestSort();

#endif
