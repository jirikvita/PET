/* File          : TSinogram.cpp
 *
 * Created       : Wed Oct 31 14:53:32 CDT 2007
 * Author        : kvita
 * Purpose       : 
 * Last modified : 
 * Comments      : For FFT in ROOT, see $ROOTSYS/cint/include/fft.c
 */

#include <string>
#include <vector>

#include "TSinogram.hpp"

// _____________________________________________
// _____________________________________________
// _____________________________________________


TSinogram::TSinogram(int binsx, double x1, double x2, int binsy, double y1, double y2) 
{
  Init(binsx, x1, x2, binsy, y1, y2, 0);
}


TSinogram::TSinogram(int binsx, double x1, double x2, int binsy, double y1, double y2, TF2 *slice) 
{
  Init(binsx, x1, x2, binsy, y1, y2, slice);
}

void TSinogram::Init(int binsx, double x1, double x2, int binsy, double y1, double y2, TF2 *slice) 
{
  _sinogram = 0;
  _reconstructed = 0;
  _slice_f2 = 0;
  _binned_slice_h2 = 0;
  _backprojection = 0;
  _sinogram_FFT1D_re = 0;
  _sinogram_FFT1D_im = 0;
  //  _sinogram_FFT2D = 0;
  _sinogram_InvFFT1D_re = 0;
  _sinogram_InvFFT1D_im = 0;
  _rand = 0;

  for (int i = 0; i < maxNbins; ++i) {
    _sinogram_re[i] = 0;
    _sinogram_im[i] = 0;
    _sinogram_Inv_re[i] = 0;
    _sinogram_Inv_im[i] = 0;
  }

  _binsx = binsx;
  _x1 = x1;
  _x2 = x2;
  _binsy = binsy;
  _y1 = y1;
  _y2 = y2;
  _slice_f2 = slice;

}
  


TSinogram::~TSinogram()
{


}

// methods

// _____________________________________________

void TSinogram::SetOriginalSlice(TF2 *slice)
{
  _slice_f2 = slice;
}

// _____________________________________________

void TSinogram::MakeTheProjection(int phibins, int sbins, double dzSF, int nsubsbins)
{

  if (!_slice_f2) {
    cout << "ERROR: no slice TF2 specified!" << endl;
    return;
  }
    
  if (_sinogram) { 
    cout << "TSinogram::MakeTheProjection: WARNING: replacing the previous sinogram..." << endl;
    delete _sinogram;
  }

  // need the diagonal of the picture:

  double SliceXmax = _slice_f2 -> GetXaxis() -> GetXmax();
  double SliceXmin = _slice_f2 -> GetXaxis() -> GetXmin();
  double SliceYmax = _slice_f2 -> GetYaxis() -> GetXmax();
  double SliceYmin = _slice_f2 -> GetYaxis() -> GetXmin();

  double Deltax = SliceXmax - SliceXmin;
  double Deltay = SliceYmax - SliceYmin;
  double diag = sqrt(  pow(Deltax, 2) + pow(Deltay, 2) );

  double mins = -diag;
  double maxs = diag;

  cout << "TSinogram::MakeTheProjection: Initialising the sinogram..." << endl;

  _sinogram = new TH2D(Form("sinogram_%s", _slice_f2 -> GetName()),
		       Form("sinogram_%s", _slice_f2 -> GetName()),
		       sbins, mins, maxs,
		       phibins, 0, TMath::Pi() );


  cout << "TSinogram::MakeTheProjection: Making the sinogram, this may take a while..." << endl;

  double ds = (maxs - mins)/sbins;  
  double dz = dzSF*ds;
  for (int i = 0; i < phibins; ++i) {
    double phi = i*TMath::Pi()/phibins;

    for (int j = 0; j < sbins; ++j) {
      for (int subj = 0; subj < nsubsbins; ++subj) {
	double s = diag - ds*(subj/(1.*nsubsbins) + j);
	// let's integrate over rotated axis z = -x sin(phi) + y cos(phi):
	double sum = 0.;
	// temporary, assuming no shift:
	double maxz = diag;
	double minz = -diag;
	int nz = (maxz - minz) / dz;

	for (int k = 0; k < nz; ++k) {
	  double z = maxz - k*dz;
	  double x1 = s*cos(phi) - z*sin(phi);
	  double y1 = s*sin(phi) + z*cos(phi);
	  sum += _slice_f2 -> Eval(x1, y1) * ds * dz;
	  /*
	    double x2 = (s+ds)*cos(phi) - (z+dz)*sin(phi);
	    double y2 = (s+ds)*sin(phi) + (z+dz)*cos(phi);
	    double xmin = (x1 < x2) ? x1 : x2;
	    double xmax = (x1 < x2) ? x2 : x1;
	    double ymin = (y1 < y2) ? y1 : y2;
	    double ymax = (y1 < y2) ? y2 : y1;
	    if (xmax < SliceXmax && xmin > SliceXmin && ymax < SliceYmax && ymin > SliceYmin) {
	    sum += _slice_f2 -> Integral(xmin, xmax, ymin, ymax, 1.e-4);
	    }
	  */
	
	} // z
	_sinogram -> SetBinContent(j+1, i+1, sum);
      } // subj
    } // j
  } // i
  

}


// _____________________________________________

void TSinogram::MakeTheBackProjection(TH2D *sinogram)
{


 if (!sinogram) {
    cout << "ERROR: sinogram is a zero pointer!" << endl;
    return;
  }
    
  if (_backprojection) { 
    cout << "TSinogram::MakeTheProjection: WARNING: replacing the previous sinogram..." << endl;
    delete _backprojection;
  }

  cout << "TSinogram::MakeTheProjection: Initialising the backprojection..." << endl;

  _backprojection = new TH2D(Form("backprojection_%s", sinogram -> GetName()),
			     Form("backprojection_%s", sinogram -> GetName()),
			     _binsx, _x1, _x2, _binsy, _y1, _y2 );


  cout << "TSinogram::MakeTheProjection: Making the backprojection, this may take a while..." << endl;
  
  for (int i = 0; i < _binsx; ++i) {
    double x = _x1 + i*(_x2 - _x1) / _binsx;
    for (int j = 0; j < _binsy; ++j) {
      double y = _y1 + j*(_y2 - _y1) / _binsy;
      double sum = 0.;
      // loop over phi's:
      for (int k = 1; k <= sinogram -> GetNbinsY(); ++k) {
	double phi = sinogram -> GetYaxis() -> GetBinCenter(k);
	//	double s = x*cos(phi) + y*sin(phi);
	double z = -x*sin(phi) + y*cos(phi);
	double delta = sinogram -> GetBinContent(sinogram -> FindBin(z, phi));
	sum += delta;
      } // phi
      sum /= sinogram -> GetNbinsX();
      _backprojection -> SetBinContent(i+1, j+1, sum);
      
    } // j
  } // i
  

}


// _____________________________________________

// supported only for power = 1 2 4 8 16 32,
// i.e. for N = 2 4 16 256 65536 2^32
int TSinogram::GetReversedBitOrder(int m, int power) 
{
  // http://www.sjbaker.org/steve/software/cute_code.html
  int n = m;
  /*
    n = ((n >>  1) & 0x55555555) | ((n <<  1) & 0xaaaaaaaa);
    n = ((n >>  2) & 0x33333333) | ((n <<  2) & 0xcccccccc);
    n = ((n >>  4) & 0x0f0f0f0f) | ((n <<  4) & 0xf0f0f0f0);
    n = ((n >>  8) & 0x00ff00ff) | ((n <<  8) & 0xff00ff00);
    n = ((n >> 16) & 0x0000ffff) | ((n << 16) & 0xffff0000);
  */

  n = ((n >>  1) & 0x55555555) | ((n <<  1) & 0xaaaaaaaa);
  if (power > 2) {
    n = ((n >>  2) & 0x33333333) | ((n <<  2) & 0xcccccccc);
    if (power > 4) {
      n = ((n >>  4) & 0x0f0f0f0f) | ((n <<  4) & 0xf0f0f0f0);
      if (power > 8) {
	n = ((n >>  8) & 0x00ff00ff) | ((n <<  8) & 0xff00ff00);
	if (power > 16) {
	  n = ((n >> 16) & 0x0000ffff) | ((n << 16) & 0xffff0000);
	}
      }
    }
  }
  // cout << "m=" << m << " n=" << n << endl;
  return n;

}

// _____________________________________________
bool TSinogram::CheckPower(int FTN)
{
  return (FTN == 1 || FTN == 2 || FTN == 8 || FTN == 16 || FTN == 32);
}

int TSinogram::CheckPowerBinning(int base, int n) 
{

  // Test to see if a number is an exact power of two:
  // bool is_power_of_two = ( (n & (n - 1)) == 0) ;

  // check that n = base^someInt

  int FTN = 0;
  int index = 1;
  while (index < n) {
    FTN++;
    index *= base;
  }
  
  cout << "Checked: n=" << n << " = " << base << "^" << FTN << endl;

  if (index == n)
    return FTN;
  else
    return -1;

}

// _____________________________________________

void TSinogram::MakeSinogramFFT1D()
{

  // go over phi bins:

  int nphi = _sinogram -> GetNbinsY();
  int nz = _sinogram -> GetNbinsX();

  // check that the binning is 2^FTN:
  int FTN = CheckPowerBinning(2, nz);
  assert(FTN > 0);
  assert(CheckPower(FTN));


  int maxFT = nz;
  int minFT = 0;

  if (_sinogram_FFT1D_re) {
    cout << "TSinogram::MakeSinogramFFT1D: WARNING: Resetting the sinogram Re 1D FFT..." << endl;
  }
  cout << "TSinogram::MakeSinogramFFT1D: Creating new sinogram Re 1D FFT ..." << endl;
  _sinogram_FFT1D_re = new TH2D(Form("FFT1D_re_%s", _sinogram -> GetName()), 
				Form("FFT1D_re_%s", _sinogram -> GetName()), nz, minFT, maxFT, 
				nphi,
				_sinogram -> GetYaxis() -> GetXmin(),
				_sinogram -> GetYaxis() -> GetXmax() );
  
  if (_sinogram_FFT1D_im) {
    cout << "TSinogram::MakeSinogramFFT1D: WARNING: Resetting the sinogram Im 1D FFT..." << endl;
  }
  cout << "TSinogram::MakeSinogramFFT1D: Creating new sinogram Im 1D FFT..." << endl;
  _sinogram_FFT1D_im = new TH2D(Form("FFT1D_im_%s", _sinogram -> GetName()), 
				Form("FFT1D_im_%s", _sinogram -> GetName()), nz, minFT, maxFT, 
				nphi,
				_sinogram -> GetYaxis() -> GetXmin(),
				_sinogram -> GetYaxis() -> GetXmax() );
  
  for (int i = 0; i < nphi; ++i) {

    if (_sinogram_re[i]) {
      cout << "Deleting the Re array..." << endl;
      delete [] _sinogram_re[i];
    }
    if (_sinogram_im[i]) {
      cout << "Deleting the Im array..." << endl;
      delete [] _sinogram_im[i];
    }

    cout << "Initialising the Re array..." << endl;
    _sinogram_re[i] = new double[nz];
    cout << "Initialising the Im array..." << endl;
    _sinogram_im[i] = new double[nz];
    
    for (int j = 0; j < nz; ++j) {
      _sinogram_re[i][j] =  _sinogram -> GetBinContent(j+1, i+1);
      _sinogram_im[i][j] = 0;
      //      cout << "N=" << FTN << " i,j=" << i << "," << j << " re: " << _sinogram_re[i][j] << " im: " << _sinogram_im[i][j] << endl;
    }
    cout << "Doing the FFT of the phi slice " << i << endl;
    FFT(FTN, _sinogram_re[i], _sinogram_im[i]);
    
    for (int j = 0; j < nz; ++j) {
      //      cout << "N=" << N << " i,j=" << i << "," << j << " re: " << _sinogram_re[i][j] << " im: " << _sinogram_im[i][j] << endl;
      _sinogram_FFT1D_re -> SetBinContent(1 + GetReversedBitOrder(j, FTN), i+1, _sinogram_re[i][j]);
      _sinogram_FFT1D_im -> SetBinContent(1 + GetReversedBitOrder(j, FTN), i+1, _sinogram_im[i][j]);
      //      _sinogram_FFT1D_re -> SetBinContent(j+1, i+1, _sinogram_re[i][j]);
      //      _sinogram_FFT1D_im -> SetBinContent(j+1, i+1, _sinogram_im[i][j]);
    }
  } // i - phi bins

  //  delete [] _sinogram_re;
  //  delete [] _sinogram_im;

}

// _____________________________________________

void TSinogram::MakeSinogramInvFFT1D()
{

  // go over phi bins:

  int nphi = _sinogram_FFT1D_re -> GetNbinsY();
  int nz = _sinogram_FFT1D_re -> GetNbinsX();

  int nz_sino = _sinogram -> GetNbinsX();

  // check that the binning is 2^FTN:
  int FTN = CheckPowerBinning(2, nz_sino);
  assert(FTN > 0);

  int maxInvFT =_sinogram -> GetXaxis() -> GetXmax();
  int minInvFT =_sinogram -> GetXaxis() -> GetXmin();
  int nInvFT = _sinogram -> GetNbinsX();

  if (_sinogram_InvFFT1D_re) {
    cout << "TSinogram::MakeSinogramFFT1D: WARNING: Resetting the sinogram Re 1D InvFFT..." << endl;
  }
  cout << "TSinogram::MakeSinogramFFT1D: Creating new sinogram Re 1D InvFFT ..." << endl;
  _sinogram_InvFFT1D_re = new TH2D(Form("InvFFT1D_re_%s", _sinogram -> GetName()), 
				   Form("InvFFT1D_re_%s", _sinogram -> GetName()), nInvFT, minInvFT, maxInvFT, 
				   nphi,
				   _sinogram -> GetYaxis() -> GetXmin(),
				   _sinogram -> GetYaxis() -> GetXmax() );
  
  if (_sinogram_InvFFT1D_im) {
    cout << "TSinogram::MakeSinogramFFT1D: WARNING: Resetting the sinogram Im 1D InvFFT..." << endl;
  }
  cout << "TSinogram::MakeSinogramFFT1D: Creating new sinogram Im 1D InvFFT..." << endl;
  _sinogram_InvFFT1D_im = new TH2D(Form("InvFFT1D_im_%s", _sinogram -> GetName()), 
				   Form("InvFFT1D_im_%s", _sinogram -> GetName()), nInvFT, minInvFT, maxInvFT, 
				   nphi,
				   _sinogram -> GetYaxis() -> GetXmin(),
				   _sinogram -> GetYaxis() -> GetXmax() );
  
  for (int i = 0; i < nphi; ++i) {

    if (_sinogram_Inv_re[i]) {
      //      cout << "Deleting the Re array..." << endl;
      delete [] _sinogram_Inv_re[i];
    }
    if (_sinogram_Inv_im[i]) {
      //      cout << "Deleting the Im array..." << endl;
      delete [] _sinogram_Inv_im[i];
    }

    //    cout << "Initialising the Re array..." << endl;
    _sinogram_Inv_re[i] = new double[nz];
    //    cout << "Initialising the Im array..." << endl;
    _sinogram_Inv_im[i] = new double[nz];
    
    for (int j = 0; j < nz; ++j) {
      _sinogram_Inv_re[i][j] = _sinogram_re[i][j];
      _sinogram_Inv_im[i][j] = _sinogram_im[i][j];
      //      cout << "N=" << FTN << " i,j=" << i << "," << j << " re: " << _sinogram_Inv_re[i][j] << " im: " << _sinogram_Inv_im[i][j] << endl;
    }
    //    cout << "Doing the FFT of the phi slice " << i << endl;
    IFFT(FTN, _sinogram_Inv_re[i], _sinogram_Inv_im[i]);
    
    for (int j = 0; j < nz; ++j) {
      //      cout << "N=" << N << " i,j=" << i << "," << j << " re: " << _sinogram_Inv_re[i][j] << " im: " << _sinogram_Inv_im[i][j] << endl;
      _sinogram_InvFFT1D_re -> SetBinContent(j+1, i+1, _sinogram_Inv_re[i][j]);
      _sinogram_InvFFT1D_im -> SetBinContent(j+1, i+1, _sinogram_Inv_im[i][j]);
    }
  } // i - phi bins

  //  delete [] _sinogram_Inv_re;
  //  delete [] _sinogram_Inv_im;

}

// _____________________________________________

TH2D* TSinogram::GetSinogramFFT1D_re()
{
  return _sinogram_FFT1D_re;
}

TH2D* TSinogram::GetSinogramFFT1D_im()
{
  return _sinogram_FFT1D_im;
}

// _____________________________________________

TH2D* TSinogram::GetSinogramInvFFT1D_re()
{
  return _sinogram_InvFFT1D_re;
}

TH2D* TSinogram::GetSinogramInvFFT1D_im()
{
  return _sinogram_InvFFT1D_im;
}


// _____________________________________________


// _____________________________________________


double TSinogram::EvalFilter(int ifilter, int channel, int maxN)
{
  

  // cool:
  if (ifilter == 0)
    return ( 2.*channel*(sin(1.5*channel*TMath::TwoPi()/maxN)) );
  if (ifilter == 1)
    return TMath::Abs( 2.*channel*(sin(1.5*channel*TMath::TwoPi()/maxN)) );
  if (ifilter == 2)
    return ( 1.*channel*(sin(1.*channel*TMath::Pi()/maxN)) );
  if (ifilter == 3)
    return ( -channel*(channel - maxN) );
  if (ifilter == 4)
    return ( maxN/2 - TMath::Abs(maxN/2 - channel ) );

  return 1.;

}

// _____________________________________________


TH1D* TSinogram::MakeFilterHisto(int ifilter, int maxN, std::string name, std::string title)
{
  
  TH1D *h = new TH1D(name.c_str(), title.c_str(), maxN, 0, maxN);
  for (int i = 0; i < maxN; ++i)
    h -> SetBinContent(i+1, EvalFilter(ifilter, i, maxN));

  return h;

}

// _____________________________________________

void TSinogram::KillSinogramChannel(int channel, double value)
{
  
  assert(_sinogram);
  int nphi = _sinogram -> GetNbinsY();
  for (int i = 0; i < nphi; ++i)
    _sinogram -> SetBinContent(channel, i + 1, value);



}

// _____________________________________________

void TSinogram::SmearSinogram(double sigma)
{

  assert(_sinogram);
  if (!_rand)
    _rand = new TRandom3();

  int nphi = _sinogram -> GetNbinsY();
  int nz = _sinogram -> GetNbinsX();
  for (int i = 0; i < nphi; ++i) {
    for (int j = 0; j < nz; ++j) {
      _sinogram -> SetBinContent(j+1, i+1, TMath::Abs(_rand -> Gaus(0., sigma)) + _sinogram -> GetBinContent(j+1, i+1));
    }
  }


}

// _____________________________________________

void TSinogram::ApplyFilter(int ifilter)
{
  
  int nphi = _sinogram -> GetNbinsY();
  int nz = _sinogram -> GetNbinsX();
  int FTN = CheckPowerBinning(2, nz);
  assert(CheckPower(FTN));

  for (int i = 0; i < nphi; ++i) {
    for (int j = 0; j < nz; ++j) {


      int nz_reversed = GetReversedBitOrder(j, FTN); 

      // make some trasformation

      double SF = EvalFilter(ifilter, nz_reversed, nz);
      _sinogram_re[i][j] = SF*_sinogram_re[i][j];
      _sinogram_im[i][j] = SF*_sinogram_im[i][j];

      //      _sinogram_re[i][j] = TMath::Abs(_sinogram_re[i][j]);
      //      _sinogram_im[i][j] = TMath::Abs(_sinogram_im[i][j]);
 
   }
  }
}

// _____________________________________________

void TSinogram::Reconstruct()
{
 // perform the deconvolution
}

// _____________________________________________

void TSinogram::DrawSinogram()
{

}

// _____________________________________________

void TSinogram::DrawReconstructed()
{

}

// _____________________________________________

TH2D* TSinogram::GetSinogram()
{
  return _sinogram;
}

TH2D* TSinogram::GetBackProjection()
{
  return _backprojection;
}

TH2D* TSinogram::GetBinnedSlice()
{
  return _binned_slice_h2;
}


TF2* TSinogram::GetOriginalSlice()
{
  return _slice_f2;
}

TH2D* TSinogram::GetReconstructedImage()
{
  return _reconstructed;

}
