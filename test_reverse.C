



int GetReversedBitOrder(int m, int power) 
{
  // http://www.sjbaker.org/steve/software/cute_code.html
  int n = m;
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
  cout << "m=" << m << " n=" << n << endl;
  return n;

}


int test_reverse()
{


  const int N = 256;
  for (int i = 0; i < N; ++i) {
    int j = GetReversedBitOrder(i, 8);
    cout << " i: " << i 
	 << " j: " << j
	 << " j%N: " << j % N
	 << endl;
  }


}
