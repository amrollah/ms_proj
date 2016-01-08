
#include <math.h>
#include <windows.h>
#include "mex.h"

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x;
  int i, i0, n;
  double m0=0., m0_, m1=0., m1_=0., score, bestscore=0.; 
  
  if (nrhs!=1 || nlhs>1 || !mxIsDouble(prhs[0])) mexErrMsgTxt("call  jmin = otsu_local_min(hist_data)");
  
  x = mxGetPr(prhs[0]);
  n = mxGetNumberOfElements(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  
  if (n<=1) {*mxGetPr(plhs[0]) = 1; return;}
  else if (n==2) {*mxGetPr(plhs[0]) = 1+(x[1]<x[0]); return; }
  
  //Otsu's method
  for (i=0; i<n; i++) {
    m0 += x[i];
    m1 += i*x[i];
  }
  m0_ = x[0];
  for (i=1;i<n-1;i++) {
    m0_ += x[i];
    if (m0_==0.) continue;
    if (m0_>=m0) break;
    m1_ += i*x[i];
    score = m1_/m0_ - (m1-m1_)/(m0-m0_);
    score = m0_*(m0-m0_)*score*score;
    if (score>bestscore) {i0 = i; bestscore=score; }
  }
  
  //find local minimum (in case of local max go towards equilibrating mass)
  if (x[i0-1]<x[i0] && (x[i0+1]>=x[i0] || m0_>0.5*(m0+x[i0]))) {
    for (i0--; i0>0 && x[i0-1]<x[i0]; i0--);
  } else if (x[i0+1]<x[i0]) {
    for (i0++; i0<n-1 && x[i0+1]<x[i0]; i0++);
  }
  
  *mxGetPr(plhs[0]) = i0+1;
}
