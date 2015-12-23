
#include <math.h>
#include <windows.h>
#include "mex.h"

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x;
  int i, i0, n;
  double m0=0., m0_=0., m1=0., m1_=0., score, bestscore=0.; 
  
  if (nrhs!=1 || nlhs>1 || !mxIsDouble(prhs[0])) mexErrMsgTxt("call  lastleft = otsu4hist(hist_data)");
  
  x = mxGetPr(prhs[0]);
  n = mxGetNumberOfElements(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  
  for (i=0; i<n; i++) {
    m0 += x[i];
    m1 += i*x[i];
  }
  for (i=0;i<n-1;i++) {
    m0_ += x[i];
    if (m0_==0.) continue;
    if (m0_>=m0) break;
    m1_ += i*x[i];
    score = m1_/m0_ - (m1-m1_)/(m0-m0_);
    score = m0_*(m0-m0_)*score*score;
    if (score>bestscore) {i0 = i; bestscore=score; }
  }
  
  *mxGetPr(plhs[0]) = i0+1;
}
