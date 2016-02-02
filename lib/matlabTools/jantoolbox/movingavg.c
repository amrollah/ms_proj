/*
function y = movingavg(x,l,r)
%   compute moving average of data x looking l steps to the left and r steps to the right
%
% 14.5.2014, Jan Poland.
*/

#include "mex.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[] )
{
  int i, n, l, r, ll, rr;
  double *x, *y;
  double nan = mxGetNaN(), v;
  
  if (nlhs>1 || nrhs!=3) mexErrMsgTxt("usage: y = movingavg(x,l,r)");
  
  if(!mxIsDouble(prhs[0]))
    mexErrMsgTxt("parameter x must be a double matrix");
  
  x = mxGetPr(prhs[0]);
  n = mxGetNumberOfElements(prhs[0]);
  
  if(!mxIsDouble(prhs[1]) || mxGetM(prhs[1])!=1  || mxGetN(prhs[1])!=1) mexErrMsgTxt("parameter l must be a scalar");
  if(!mxIsDouble(prhs[2]) || mxGetM(prhs[2])!=1  || mxGetN(prhs[2])!=1) mexErrMsgTxt("parameter r must be a scalar");
  
  l = (int)mxGetScalar(prhs[1]);
  r = (int)mxGetScalar(prhs[2]);
  
  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),mxREAL);
  y = mxGetPr(plhs[0]);
  
  ll = max(0,-l);
  rr = min(n-1,r);
  for (v=0.0,i=ll;i<=rr;i++) v+=x[i];
  
  for (i=0;;) {
    if (rr-ll<0) *(y++)=nan; else *(y++)=v/(rr-ll+1);
    if (++i==n) break;
    if (i-1-l>=0) {ll = i-l; if (i-1-l<n) v-=x[i-1-l];}
    if (i+r<=n-1) {rr = i+r; if (i+r>=0) v+=x[i+r];}
  }
}
