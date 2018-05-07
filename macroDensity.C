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


#ifndef MACRODENSITY_C
#define MACRODENSITY_C

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

const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; 

//______________________________________________________________________________
//Double_t TStat::Mean(Int_t n, const Double_t *arr)
Double_t Mean(Int_t n, const Double_t *arr)
{
   // Calculate arithmetic mean
   if(kCSa) cout << "------TStat::Mean------" << endl;

   if (n == 1) return arr[0];

   Double_t mean = 0.0;
   for (Int_t i=0; i<n; i++) mean += arr[i];

   return mean/n;
}//Mean

//______________________________________________________________________________
//Double_t TStat::Var(Int_t n, const Double_t *arr, const Double_t mean)
Double_t Var(Int_t n, const Double_t *arr, const Double_t mean)
{
   // Calculate variance for arithmetic mean
   if(kCSa) cout << "------TStat::Var------" << endl;

   if (n == 1) return 0;

   Double_t var = 0.0;
   Double_t tmp = 0.0;
   for (Int_t i=0; i<n; i++) {
      tmp = arr[i] - mean;
      var += tmp * tmp;
   }//for_i

   return var / (n - 1);
}//Var

//______________________________________________________________________________
void MassDist(Int_t nx, Double_t *x, Double_t *w, Double_t xmin, Double_t xmax,
              Int_t ny, Double_t *y)
{
   // Discretize data x with xmin and xmax minimum and maximum value of x
   // Array w of length nx is weight for x
   // Output array y contains discretation scheme of data
   // see Applied Statistics R50 and Applied Statistics 176
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCS) cout << "------MassDist------" << endl;
  
   Int_t    ixmin = 0;
   Int_t    ixmax = ny - 2;
   Double_t mass  = 0.0;
   Double_t delta = (xmax - xmin) / (ny - 1);
  
   for (Int_t i=0; i<ny; i++) y[i] = 0.0;

   for (Int_t i=0; i<nx; i++) mass += w[i];
   mass = 1.0/mass;

   for(Int_t i=0; i<nx; i++) {
      if (TMath::Finite(x[i])) {
         Double_t pos = (x[i] - xmin) / delta;
         Int_t    ix  = (Int_t)TMath::Floor(pos);
         Double_t fx  = pos - ix;

         if ((ixmin <= ix) && (ix <= ixmax)) {
            y[ix]     += w[i]*(1 - fx);
            y[ix + 1] += w[i]*fx;
         } else if (ix == -1) {
            y[0] += w[i]*fx;
         } else if (ix == ixmax + 1) {
            y[ix] += w[i]*(1 - fx);
         }//if
      }//if
   }//for_i
  
   for (Int_t i=0; i<ny; i++) y[i] *= mass;
}//MassDist

//______________________________________________________________________________
void TwiddleFactor4FFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im)
{
   // Twiddle factor for FFT
   // n     - length of data series
   // tf_re - on output contains real part of twiddle factor
   // tf_im - on output contains imaginary part of twiddle factor
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TwiddleFactor4FFT------" << endl;

   if (i ==0 ) {
      tf_re = 1;
      tf_im = 0;
   } else {
      tf_re =  TMath::Cos(2*TMath::Pi()*(Double_t)i / (Double_t)n);  
      tf_im = -TMath::Sin(2*TMath::Pi()*(Double_t)i / (Double_t)n); 
   }//if
}//TwiddleFactor4FFT

//______________________________________________________________________________
void TwiddleFactor4IFFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im)
{
   // Twiddle factor for Inverse FFT
   // n     - length of data series
   // tf_re - on output contains real part of twiddle factor
   // tf_im - on output contains imaginary part of twiddle factor
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TwiddleFactor4IFFT------" << endl;

   if (i == 0) {
      tf_re = 1;
      tf_im = 0;
   } else {
      tf_re = TMath::Cos(2*TMath::Pi()*(Double_t)i / (Double_t)n);  
      tf_im = TMath::Sin(2*TMath::Pi()*(Double_t)i / (Double_t)n); 
   }//if
}//TwiddleFactor4IFFT

//______________________________________________________________________________
void FFT(Int_t n, Double_t *f_re, Double_t *f_im)
{
   // Fast Fourier Transform in space, result is in reverse bit order.
   // compute FFT using Decimation In Frequency of a data sequence of length 2^n
   // f_re - real component of data series
   // f_im - imaginary component of data series
   // n    -  where 2^n is length of data series
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCS) cout << "------FFT------" << endl;
  
   Int_t    baseE, baseO, blocks, points, points2;
   Double_t even_re, even_im, odd_re, odd_im, tf_re, tf_im;

   blocks = 1;
   points = 1 << n;

   for (Int_t i=0; i<n; i++) {
      points2 = points >> 1;
      baseE = 0;

      /*
	cout
	<< " points=" << points
	<< " points2=" << points2
	<< endl;
      */

      //      cout << "i=" << i << endl;

      for (Int_t j=0; j<blocks; j++) {
         baseO = baseE + points2;

	 /*
	   cout << "j=" << j 
	   << " baseE=" << baseE
	   << " baseO=" << baseO
	   << endl;
	 */

         for (Int_t k=0; k<points2; k++) {

	   /*
	     cout << "k=" << k 
	     << "baseE+k=" << baseE + k
	     << "baseO+k=" << baseO + k
	     << endl;
	   */

            even_re = f_re[baseE + k] + f_re[baseO + k]; 
            even_im = f_im[baseE + k] + f_im[baseO + k];

	    /*
	    cout
	    << " even_re=" << even_re
	    << " even_im=" << even_im
	    << endl;
	    cout
	      << " points=" << points
	      << " k=" << k
	      << " tf_re=" << tf_re
	      << " tf_im=" << tf_im
	      << endl;
	   */
	    
	    TwiddleFactor4FFT(points, k, tf_re, tf_im);

            odd_re = (f_re[baseE + k] - f_re[baseO + k])*tf_re
                   - (f_im[baseE + k] - f_im[baseO + k])*tf_im;
            odd_im = (f_re[baseE + k] - f_re[baseO + k])*tf_im
                   + (f_im[baseE + k] - f_im[baseO + k])*tf_re; 

	    /*
	      cout
	      << " odd_re=" << odd_re
	      << " odd_im=" << odd_im
	      << endl;
	    */

	    f_re[baseE + k] = even_re;
            f_im[baseE + k] = even_im;
            f_re[baseO + k] = odd_re;
            f_im[baseO + k] = odd_im;
         }//for_k

         baseE = baseE + points;
      }//for_j 
                   
      blocks = blocks << 1; 
      points = points >> 1;
   }//for_i
}//FFT

//______________________________________________________________________________
void IFFT(Int_t n, Double_t *f_re, Double_t *f_im)
{
   // Inverse Fast Fourier Transform in space, where input is in reverse bit
   // order and output is in normal order.
   // compute IFFT using Decimation In Time of a data sequence of length 2^n
   // f_re - real component of data series
   // f_im - imaginary component of data series
   // n    -  where 2^n is length of data series
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCS) cout << "------IFFT------" << endl;

   Int_t    blocks, points, points2, baseB, baseT;
   Double_t top_re, top_im, bot_re, bot_im, tf_re, tf_im;

   blocks = 1 << (n-1);
   points = 2;  
   for (Int_t i=0; i<n; i++) {
      points2 = points >> 1;
      baseT = 0;

      for (Int_t j=0; j<blocks; j++) {
         baseB = baseT + points2;

         for (Int_t k=0; k<points2; k++) {
            top_re = f_re[baseT + k];
            top_im = f_im[baseT + k];
	
            TwiddleFactor4IFFT(points, k, tf_re, tf_im);

            bot_re = f_re[baseB + k]*tf_re - f_im[baseB + k]*tf_im;
            bot_im = f_re[baseB + k]*tf_im + f_im[baseB + k]*tf_re;

            f_re[baseT + k] = top_re + bot_re;
            f_im[baseT + k] = top_im + bot_im;
            f_re[baseB + k] = top_re - bot_re;   
            f_im[baseB + k] = top_im - bot_im; 
         }//for_k  
  
         baseT = baseT + points; 
      }//for_j

      blocks = blocks >> 1;
      points = points << 1;
   }//for_i
}//IFFT

//______________________________________________________________________________
Int_t FFTDensityConvolve(Int_t n, Double_t *x_re, Double_t *y_re)
{
   // Fast Fourier Transform density convolve
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCS) cout << "------FFTDensityConvolve------" << endl;

   Int_t err = 0;

   // to stop rounding problems 
   Int_t nlog2 = (Int_t)(TMath::Log((Double_t)n) / TMath::Log(2.0) + 0.5);

// Init local arrays
   Double_t *x_im = 0;
   Double_t *y_im = 0;
   Double_t *c_re = 0;
   Double_t *c_im = 0;

// Create local arrays
   if (!(x_im = new Double_t[n])) {err = 1; goto cleanup;}
   if (!(y_im = new Double_t[n])) {err = 1; goto cleanup;}
   if (!(c_re = new Double_t[n])) {err = 1; goto cleanup;}
   if (!(c_im = new Double_t[n])) {err = 1; goto cleanup;}

   FFT(nlog2, y_re, y_im);
   FFT(nlog2, x_re, x_im);
  
   for (Int_t i=0; i<n; i++){
      c_re[i] = y_re[i]*x_re[i] + y_im[i]*x_im[i];
      c_im[i] = y_im[i]*x_re[i] - y_re[i]*x_im[i];
   }//for_i
  
   IFFT(nlog2, c_re, c_im);

   for (Int_t i=0; i<n; i++){
      x_re[i] = c_re[i];
   }//for_i

cleanup:
   delete [] c_im;
   delete [] c_re;
   delete [] y_im;
   delete [] x_im;

   return err;
}//FFTDensityConvolve

//______________________________________________________________________________
void Kernelize(Int_t n, Double_t *x, Double_t bw, const char *kernel)
{
   // Kernelize arry x of length n with bandwidth bw
   // kernel is the type of kernel to use
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ and expanded for XPS by Christian Stratowa 
   if(kCS) cout << "------Kernelize------" << endl;
  
   Double_t a, ax;

   if (strcmp(kernel, "gaussian") == 0) {
   // Gaussian Kernel
      for (Int_t i=0; i<n; i++) {
         x[i] = TMath::Gaus(x[i], 0, bw, kTRUE);
      }//for_i
   } else if (strcmp(kernel, "epanechnikov") == 0) {
   // Epanechnikov Kernel
      a = bw * TMath::Sqrt(5.0);
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = 3.0/(4.0*a)*(1.0 - (ax/a)*(ax/a));
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "rectangular") == 0) {
   // Rectangular Kernel
      a = bw * TMath::Sqrt(3.0);
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = 0.5 / a;
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "triangular") == 0) {
   // Triangular Kernel
      a = bw * TMath::Sqrt(6.0);
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = (1.0 - ax/a) / a;
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "biweight") == 0) {
   // Biweight Kernel
      a = bw * TMath::Sqrt(7.0);
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = 15.0/16.0 * (1.0 - (ax/a)*(ax/a)) * (1.0 - (ax/a)*(ax/a)) / a;
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "cosine") == 0) {
   // Cosine Kernel
      a = bw / TMath::Sqrt(1.0/3.0 - 2.0/(TMath::Pi()*TMath::Pi()));
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = (1.0 + TMath::Cos(TMath::Pi()*x[i]/a)) / (2*a);
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "optcosine") == 0) {
   // Optcosine Kernel
      a = bw / TMath::Sqrt(1.0 - 8.0/(TMath::Pi()*TMath::Pi()));
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = (TMath::Pi()/4.0)*TMath::Cos(TMath::Pi()*x[i]/(2*a)) / a;
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   }//if
}//Kernelize

//______________________________________________________________________________
Double_t Bandwidth(Int_t n, Double_t *x, Double_t iqr)
{
   // Compute kernel bandwidth for array x of size n
   // iqr is the interquartile range of x
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCS) cout << "------Bandwidth------" << endl;

   Double_t hi, lo;
  
//   hi = TMath::RMS(n, x); // (x-m)/n not (x-m)/(n-1) =>is not stdev!!
   hi = TMath::Sqrt(Var(n, x, Mean(n, x)));
cout << "hi = " << hi << endl; 
  
   if (hi > iqr) lo = iqr/1.34;
   else          lo = hi;

   if (lo == 0) {
      if      (hi != 0)               lo = hi;
      else if (TMath::Abs(x[1]) != 0) lo = TMath::Abs(x[1]);
      else                            lo = 1.0;
   }//if

   return (0.9*lo*TMath::Power((Double_t)n, -0.2));
}//Bandwidth 

//______________________________________________________________________________
void LinearInterpolate(Double_t *xin, Double_t *yin,
                       Int_t nout, Double_t *xout, Double_t *yout)
{
   // Given xin and yin, interpolate linearly at xout and put results in yout
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCS) cout << "------LinearInterpolate------" << endl;

   Int_t i, j, ij;

   for (Int_t k=0 ; k<nout; k++) {
      i = 0;
      j = nout - 1;
    
      if(xout[k] < xin[i]) {yout[k] = yin[0];      continue;}
      if(xout[k] > xin[j]) {yout[k] = yin[nout-1]; continue;}
 
      // find the correct interval by bisection
      while (i < j - 1) {  // xin[i] <= xout[k] <= xin[j]
         ij = (i + j)/2;   // i+1 <= ij <= j-1
         if (xout[k] < xin[ij]) j = ij;
         else                   i = ij; // still i < j
      }//while
  
      if (xout[k] == xin[j]) {yout[k] = yin[j]; continue;}
      if (xout[k] == xin[i]) {yout[k] = yin[i]; continue;}
  
      yout[k] = yin[i] + (yin[j] - yin[i])*((xout[k] - xin[i])/(xin[j] - xin[i]));
   }//for_k
}//LinearInterpolate

//______________________________________________________________________________
Int_t Density(Int_t n, Double_t *x, Double_t *w, Int_t nout,  Double_t *xout,
              Double_t *yout, const char *kernel = "epanechnikov")
{
   // Kernel density of array x with weight w of size n
   // yout is the array of density values of size nout
   // xout is the array of the corresponding x coordinates
   // Note: size nout of output should be a power of two, preferably 512 or above
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCS) cout << "------Density------" << endl;

   Int_t err   = 0;
   Int_t nout2 = 2*nout;
   Double_t lo, hi, iqr, bw, from, to;

// Init local arrays
   Int_t    *indx = 0;
   Double_t *sort = 0;
   Double_t *xin  = 0;
   Double_t *yin  = 0;
   Double_t *yord = 0;

// Create local arrays
   if (!(indx = new Int_t[n]))        {err = 1; goto cleanup;}
   if (!(sort = new Double_t[n]))     {err = 1; goto cleanup;}
   if (!(xin  = new Double_t[nout]))  {err = 1; goto cleanup;}
   if (!(yin  = new Double_t[nout2])) {err = 1; goto cleanup;}
   if (!(yord = new Double_t[nout2])) {err = 1; goto cleanup;}

// Sort x to get lo, hi and interquartile range
//??   TMath::Sort(n, x, indx);
   TMath::Sort(n, x, indx, kFALSE);
   for (Int_t i=0; i<n; i++) sort[i] = x[indx[i]];
cout << "sort[0] = " << sort[0] << endl; 
cout << "sort[n-1] = " << sort[n-1] << endl; 
  
   lo  = sort[0];
   hi  = sort[n-1];
   iqr = sort[(Int_t)(0.75*n + 0.5)] - sort[(Int_t)(0.25*n + 0.5)];
cout << "lo = " << lo << endl; 
cout << "hi = " << hi << endl; 
cout << "iqr = " << iqr << endl; 

// Bandwidth  
   bw = Bandwidth(n, x, iqr);
cout << "bw = " << bw << endl; 
  
   lo = lo - 7*bw;
   hi = hi + 7*bw;

   for (Int_t i=0; i<=nout; i++) {
      yin[i] = (Double_t)i/(Double_t)(nout2 - 1)*2*(hi - lo);
   }//for_i

   for (Int_t i=nout+1; i<nout2; i++) {
      yin[i] = -yin[nout2 - i];
   }//for_i

   Kernelize(nout2, yin, bw, kernel);

   MassDist(n, x, w, lo, hi, nout, yord);

   err = FFTDensityConvolve(nout2, yin, yord);
   if (err != 0) goto cleanup;

   // corrections to get correct output range 
   to   = hi - 4*bw;
   from = lo + 4*bw;

   for (Int_t i=0; i<nout; i++) {
      xin[i]  = (Double_t)i / (Double_t)(nout - 1)*(hi - lo)   + lo;
      xout[i] = (Double_t)i / (Double_t)(nout - 1)*(to - from) + from;
   }//for_i

   for (Int_t i=0; i<nout; i++) {
      yin[i] = yin[i]/nout2;
   }//for_i

   LinearInterpolate(xin, yin, nout, xout, yout);

cleanup:
   delete [] xin;
   delete [] yord;
   delete [] yin;
   delete [] sort;
   delete [] indx;

   return err;
}//Density

//______________________________________________________________________________
void TestDensity(const char *kernel = "epanechnikov")
{
   if(kCS) cout << "------TestDensity------" << endl;

   Int_t npts = 512;
   Int_t nrow = 20;
   Double_t x[] = {-20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20};

   Double_t *weights = new Double_t[nrow];
   Double_t *dens_x  = new Double_t[npts];
   Double_t *dens_y  = new Double_t[npts];

   for (Int_t j=0; j<nrow; j++) {
      weights[j] = 1.0;
   }
  
   Density(nrow, x, weights, npts, dens_x, dens_y, kernel);
cout << "npts = " << npts << endl; 

   TCanvas *fCanvas = new TCanvas("test", "test", 10, 10, 400, 400);
   TGraph  *graph   = new TGraph(npts, dens_x, dens_y);
   graph->Draw("ACP");
   fCanvas->Update();

//   delete [] dens_y;
//   delete [] dens_x;
   delete [] weights;
}//TestDensity

//______________________________________________________________________________
void TestSort()
{
   Int_t n = 20;
   Double_t x[] = {-20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20};

   Int_t    *indx = new Int_t[n];
   Double_t *sort = new Double_t[n];

   TMath::Sort(n, x, indx, kFALSE);
   for (Int_t i=0; i<n; i++) sort[i] = x[indx[i]];
cout << "sort[0] = " << sort[0] << endl; 
cout << "sort[n-1] = " << sort[n-1] << endl; 
  

   delete [] sort;
   delete [] indx;
}//TestDensity

#endif
