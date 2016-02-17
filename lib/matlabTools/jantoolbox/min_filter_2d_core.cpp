
#include <windows.h>
#include <math.h>
#include <stdint.h>
#include "mex.h"

using namespace std;

int m,n,circ_mask,r,zero_is_nan;
double min_def_frac, nan;
int *hist, *ry;

template <typename T> void core(T *x, T *y, int nhist) {
  int i,j,a,c,s,ii,jj,jjmin,jjmax,iimax;
  
  for(j=0;j<n;j++) {
    for (i=0;i<nhist;i++) hist[i]=0;
    jjmin = max(0,j-r); jjmax = min(n,j+1+r);
    for (a=0,jj=jjmin;jj<jjmax;jj++) {
      if (circ_mask) iimax=min(m,ry[jj-j+r]); else iimax=min(m,r);
      for (ii=0;ii<iimax;ii++) {a++; hist[x[jj*m+ii]]++;}
    }
    if (circ_mask) {
      for (i=0;i<m;i++,y++) {
        for (jj=jjmin;jj<jjmax;jj++)
          if ((ii=i-1-ry[jj-j+r])>=0) {a--; hist[x[jj*m+ii]]--;}
        for (jj=jjmin;jj<jjmax;jj++)
          if ((ii=i+ry[jj-j+r])<m) {a++; hist[x[jj*m+ii]]++;}
          
        if (zero_is_nan) {c=a-hist[0];ii=1;} else {c=a;ii=0;}
        if ((double)c/a<min_def_frac) *y=0; else { //locate minimum
          for(s=hist[ii];s<1;s+=hist[++ii]);
          *y=(T)(ii);
        }
      }
    } else {
      for (i=0;i<m;i++,y++) {
        if ((ii=i-1-r)>=0)
          for (jj=jjmin;jj<jjmax;jj++) {a--; hist[x[jj*m+ii]]--;}
        if ((ii=i+r)<m)
          for (jj=jjmin;jj<jjmax;jj++) {a++; hist[x[jj*m+ii]]++;}

        if (zero_is_nan) {c=a-hist[0];ii=1;} else {c=a;ii=0;}
        if ((double)c/a<min_def_frac) *y=0; else { //locate minimum
          for(s=hist[ii];s<1;s+=hist[++ii]);
          *y=(T)(ii);
        }
        
      }
    }
  }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double v;
  int i;
  mxClassID classid;
  
  if (nrhs!=5 || nlhs>1) mexErrMsgTxt("call  y = min_filter_2d_core(x,r,zero_is_nan,circ_mask,max_nan_frac)");
  for (i=1;i<nrhs;i++) if (!mxIsDouble(prhs[i])) mexErrMsgTxt("r, zero_is_nan, circ_mask, and max_nan_frac must be double scalars");
  for (i=1;i<nrhs;i++) if (mxGetM(prhs[i])!=1 || mxGetN(prhs[i])!=1) mexErrMsgTxt("r, zero_is_nan, circ_mask, and max_nan_frac must be double scalars");
  
  classid = mxGetClassID(prhs[0]);
  m = (int)mxGetM(prhs[0]);
  n = (int)mxGetN(prhs[0]);
  if (classid!=mxUINT8_CLASS && classid!=mxUINT16_CLASS) mexErrMsgTxt("x must be of class uint8 or uint16");
  
  zero_is_nan = (int)mxGetScalar(prhs[2]);
  circ_mask = (int)mxGetScalar(prhs[3]);
  min_def_frac = 1.-mxGetScalar(prhs[4]);
  
  v = mxGetScalar(prhs[1]); r = (int)floor(v);
  
  if (circ_mask) {
    v = v*v;
    ry = (int*)mxMalloc((2*r+1)*sizeof(int));
    ry[r]=r;
    for (i=1;i<=r;i++) {ry[r+i] = (int)floor(sqrt(v-i*i)); ry[r-i]=ry[r+i];}
  }
  
  nan = mxGetNaN();
  
  plhs[0] = mxCreateNumericMatrix(m, n, classid, mxREAL);

  if (classid==mxUINT8_CLASS) {
    hist = (int*)mxMalloc(256*sizeof(int));
    core((uint8_t*)mxGetData(prhs[0]),(uint8_t*)mxGetData(plhs[0]),256);
  } else {
    hist = (int*)mxMalloc(65536*sizeof(int));
    core((uint16_t*)mxGetData(prhs[0]),(uint16_t*)mxGetData(plhs[0]),65536);
  }
  
  mxFree(hist);
  if (circ_mask) mxFree(ry);
}
