
#include <math.h>
#include <windows.h>
#include "mex.h"

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *y, *sep;
  int *hist;
  double resolution, score, bestscore, ymin, ymax, sum, sum1;
  int n, bin0, nbins, i, w, j0, j1;
  
  if (nrhs!=2 || nlhs>1) mexErrMsgTxt("call ysep = otsu(y, resolution)");
  
  if (!mxIsDouble(prhs[0])) mexErrMsgTxt("y must contain doubles");
  if (!mxIsDouble(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1) 
    mexErrMsgTxt("resolution must be scalar");
  
  y = mxGetPr(prhs[0]);
  n = (int)mxGetNumberOfElements(prhs[0]);
  resolution = mxGetScalar(prhs[1]);
  
  ymin = y[0]; ymax = y[0];
  for (i=1;i<n;i++) {
    if (y[i]<ymin) ymin=y[i];
    if (y[i]>ymax) ymax=y[i];
  }
  ymin = floor(.5+ymin/resolution);
  ymax = floor(.5+ymax/resolution);
  bin0 = (int)ymin;
  nbins = 1+(int)ymax-bin0;
  
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  sep = mxGetPr(plhs[0]);
  
  //compute data histogram
  hist = (int*)mxMalloc(nbins*sizeof(int));
  for (i=0;i<nbins;i++) hist[i]=0;
  for (i=0;i<n;i++) hist[(int)floor(.5+y[i]/resolution)-bin0]++;
  
  //Otsu's method for separating two peaks
  for (i=1, sum=0.; i<nbins; ++i) sum += i * hist[i];
  bestscore = 0.; sum1 = 0.; w = 0;
  for (i=0;i<nbins;i++) {
    w += hist[i];
    if (w==0) continue;
    if (w==n) break;
    sum1 += i*hist[i];
    score = sum1/w - (sum-sum1)/(n-w);
    score = score*score*w*(n-w);
    if (score>bestscore) {j0 = i; j1 = i; bestscore=score; } 
    else if (score==bestscore) j1 = i;
  }
  *sep = (0.5*(j0+j1+1)+bin0)*resolution;

  mxFree(hist);
}
