/*
function d = maxabsdist(x,y)
%   calculate the max abs distance matrix x-y (or x-x)
%
% 28.1.2012, Jan Poland.
*/

#include "mex.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[] )
{
  int i, j, k, n, m, d;
  double *x, *y, *z, s, v;
  double eps = mxGetEps(), inf = mxGetInf();
  
  if (nrhs>2) mexErrMsgTxt("Too many input arguments");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments");
  
  if(nrhs<1 || !mxIsDouble(prhs[0]) || mxIsEmpty(prhs[0]))
    mexErrMsgTxt("Parameter x must be a double matrix");
  
  x = mxGetPr(prhs[0]);  
  d = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  
  if (nrhs<2) {
    y = x;
    m = n;
  } else {
    
    if(!mxIsDouble(prhs[1]) || mxGetM(prhs[1])!=d)
      mexErrMsgTxt("Parameter y must have correct dimension");
        
    y = mxGetPr(prhs[1]);  
    m = mxGetN(prhs[1]);
  }
  
  plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL);
  
  z = mxGetPr(plhs[0]);
  
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) {
      s = 0.0;
      for (k=0;k<d;k++) {
        v = x[j*d+k]-y[i*d+k];
        if (v>s) s=v;
        else if (-v>s) s = -v;
      }
      *(z++) = s;
    }
  }
}
