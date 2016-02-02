
#include <windows.h>
#include <math.h>
#include "mex.h"

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x, *y, *x1;
  double v, s, min_def_frac, nan;
  int ksize,m,n,i,j,c,a,ii,jj,jjmin,jjmax;
  
  if (nrhs<2 || nrhs>3 || nlhs>1) mexErrMsgTxt("call  y = moving_avg_2d(x,ksize,[max_nan_frac = 0.5])");
  for (i=0;i<nrhs;i++) if (!mxIsDouble(prhs[i])) mexErrMsgTxt("all arguments must be double");
  for (i=1;i<nrhs;i++) if (mxGetM(prhs[i])!=1 || mxGetN(prhs[i])!=1) mexErrMsgTxt("ksize and max_nan_frac must be scalar");
  
  x = mxGetPr(prhs[0]);
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  
  v = mxGetScalar(prhs[1]);
  if (v<1. || (v-1.)/2.-floor((v-1.)/2.)!=0.) mexErrMsgTxt("ksize must be an odd positive integer");
  ksize = (int)((v-1.)/2.);
  
  if (nrhs>=3) min_def_frac = 1.-mxGetScalar(prhs[2]); else min_def_frac = 0.5;
  
  nan = mxGetNaN();
  
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
  y = mxGetPr(plhs[0]);
  
  if (m>=n) { //go column-wise
    for(j=0;j<n;j++) {
      jjmin = max(0,j-ksize); jjmax = min(n,j+1+ksize);
      for (s=0.,c=0,a=0,jj=jjmin;jj<jjmax;jj++)
        for (ii=0;ii<min(m,ksize);ii++) {
            v = x[jj*m+ii]; a++;
            if (!mxIsNaN(v)) {c++; s+=v;}
          }
      for (i=0;i<m;i++,y++) {
        if ((ii=i-1-ksize)>=0)
          for (jj=jjmin;jj<jjmax;jj++) {
            v = x[jj*m+ii]; a--;
            if (!mxIsNaN(v)) {c--; s-=v;}
          }
        if ((ii=i+ksize)<m)
          for (jj=jjmin;jj<jjmax;jj++) {
            v = x[jj*m+ii]; a++;
            if (!mxIsNaN(v)) {c++; s+=v;}
          }
        if ((double)c/a>min_def_frac) *y=s/c; else *y=nan;
      }
    }
  } else { //go row-wise
    for(j=0;j<m;j++) {
      jjmin = max(0,j-ksize); jjmax = min(m,j+1+ksize);
      for (s=0.,c=0,a=0,jj=jjmin;jj<jjmax;jj++)
        for (ii=0;ii<min(n,ksize);ii++) {
            v = x[ii*m+jj]; a++;
            if (!mxIsNaN(v)) {c++; s+=v;}
          }
      for (i=0;i<n;i++) {
        if ((ii=i-1-ksize)>=0)
          for (jj=jjmin;jj<jjmax;jj++) {
            v = x[ii*m+jj]; a--;
            if (!mxIsNaN(v)) {c--; s-=v;}
          }
        if ((ii=i+ksize)<n)
          for (jj=jjmin;jj<jjmax;jj++) {
            v = x[ii*m+jj]; a++;
            if (!mxIsNaN(v)) {c++; s+=v;}
          }
        if ((double)c/a>min_def_frac) y[i*m+j]=s/c; else y[i*m+j]=nan;
      }
    }
  }
}
