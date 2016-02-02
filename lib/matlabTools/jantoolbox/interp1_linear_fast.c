/*
function z = interp1_linear_fast(xa,ya,xx)
% yy = interp1_fast(xa,ya,xx)
%  Fast linear interpolation in 1d. Attention: xa and xx must be sorted increasingly.
*/

#include <math.h>
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[] )
{
  int i, j=0, n, m;
  double *xa, *ya, *xx, *yy;
  double d;
  
  if (nrhs!=3) mexErrMsgTxt("call yy = interp1_fast(xa,ya,xx)");
  if (nlhs>1) mexErrMsgTxt("call yy = interp1_fast(xa,ya,xx)");
    
  if(!mxIsDouble(prhs[0])) mexErrMsgTxt("xa must contain doubles"); 
  n = mxGetNumberOfElements(prhs[0]);
  xa = mxGetPr(prhs[0]);  
  
  if(!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1])!=n) mexErrMsgTxt("ya must contain doubles and match xa in size"); 
  n = mxGetNumberOfElements(prhs[0]);
  ya = mxGetPr(prhs[1]);  
  
  if(!mxIsDouble(prhs[2])) mexErrMsgTxt("xx must contain doubles"); 
  m = mxGetNumberOfElements(prhs[2]);
  xx = mxGetPr(prhs[2]);  
  
  //no checks in order to be really fast
  //for (i=0;i<n-1;i++) if (xa[i+1]<xa[i]) mexErrMsgTxt("xa must be sorted increasingly"); 
  //for (i=0;i<m-1;i++) if (xx[i+1]<xx[i]) mexErrMsgTxt("xx must be sorted increasingly"); 
  
  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[2]),mxGetN(prhs[2]),mxREAL);
  yy = mxGetPr(plhs[0]);
  
  for (i=0; i<m; i++)
  {
    while (j<n && xa[j]<xx[i]) j++;
    if (j==0) yy[i] = ya[j];
    else if (j==n) yy[i] = ya[n-1];
    else {
      d = xa[j]-xa[j-1];
      if (d<1.e-15) yy[i]=xa[j];
      else {
        d = (xx[i]-xa[j-1])/d;
        yy[i]=(1.0-d)*ya[j-1]+d*ya[j];
      }
    }
  }
}
