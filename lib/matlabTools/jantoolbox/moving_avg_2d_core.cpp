
#include <windows.h>
#include <math.h>
#include "mex.h"

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x, *y;
  int *ry;
  double v, s, min_def_frac, nan;
  int m,n,i,j,c,a,ii,jj,jjmin,jjmax,iimax,circ_mask,r;
  
  if (nrhs!=4 || nlhs>1) mexErrMsgTxt("call  y = moving_avg_2d_core(x,r,circ_mask,max_nan_frac)");
  for (i=0;i<nrhs;i++) if (!mxIsDouble(prhs[i])) mexErrMsgTxt("all arguments must be double");
  for (i=1;i<nrhs;i++) if (mxGetM(prhs[i])!=1 || mxGetN(prhs[i])!=1) mexErrMsgTxt("r, circ_mask, and max_nan_frac must be scalars");
  
  x = mxGetPr(prhs[0]);
  m = (int)mxGetM(prhs[0]);
  n = (int)mxGetN(prhs[0]);
  
  circ_mask = (int)mxGetScalar(prhs[2]);
  min_def_frac = 1.-mxGetScalar(prhs[3]);
  
  v = mxGetScalar(prhs[1]); r = (int)floor(v);
  
  if (circ_mask) {
    v = v*v;
    ry = (int*)mxMalloc((2*r+1)*sizeof(int));
    ry[r]=r;
    for (i=1;i<=r;i++) {ry[r+i] = (int)floor(sqrt(v-i*i)); ry[r-i]=ry[r+i];}
  }
  
  nan = mxGetNaN();
  
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
  y = mxGetPr(plhs[0]);
  
  for(j=0;j<n;j++) {
    jjmin = max(0,j-r); jjmax = min(n,j+1+r);
    for (s=0.,c=0,a=0,jj=jjmin;jj<jjmax;jj++) {
      if (circ_mask) iimax=min(m,ry[jj-j+r]); else iimax=min(m,r);
      for (ii=0;ii<iimax;ii++) {
        v = x[jj*m+ii]; a++;
        if (!mxIsNaN(v)) {c++; s+=v;}
      }
    }
    if (circ_mask) {
      for (i=0;i<m;i++,y++) {
        for (jj=jjmin;jj<jjmax;jj++)
          if ((ii=i-1-ry[jj-j+r])>=0) {
            v = x[jj*m+ii]; a--;
            if (!mxIsNaN(v)) {c--; s-=v;}
          }
        for (jj=jjmin;jj<jjmax;jj++)
          if ((ii=i+ry[jj-j+r])<m) {
            v = x[jj*m+ii]; a++;
            if (!mxIsNaN(v)) {c++; s+=v;}
          }
        if ((double)c/a>=min_def_frac) *y=s/c; else *y=nan;
      }
    } else {
      for (i=0;i<m;i++,y++) {
        if ((ii=i-1-r)>=0)
          for (jj=jjmin;jj<jjmax;jj++) {
            v = x[jj*m+ii]; a--;
            if (!mxIsNaN(v)) {c--; s-=v;}
          }
        if ((ii=i+r)<m)
          for (jj=jjmin;jj<jjmax;jj++) {
            v = x[jj*m+ii]; a++;
            if (!mxIsNaN(v)) {c++; s+=v;}
          }
        if ((double)c/a>=min_def_frac) *y=s/c; else *y=nan;
      }
    }
  }

  if (circ_mask) mxFree(ry);
}
