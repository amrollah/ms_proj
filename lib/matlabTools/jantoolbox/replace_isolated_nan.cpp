
#include <windows.h>
#include <math.h>
#include "mex.h"

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x, *y, *x1;
  double v, v1;
  int m,n,i,j,c,a,ii,jj;
  
  if (nrhs!=1 || nlhs>1 || !mxIsDouble(prhs[0])) mexErrMsgTxt("call  y = replace_isolated_nan(x)");
  
  x = mxGetPr(prhs[0]);
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
  y = mxGetPr(plhs[0]);
  
  for(x1=x,j=0;j<n;j++) 
    for (i=0;i<m;i++,x1++,y++) {
      if (!mxIsNaN(*x1)) *y=*x1; else {
        for (c=0,a=0,v=0.,jj=max(0,j-1);jj<min(n,j+2);jj++)
          for (ii=max(0,i-1);ii<min(m,i+2);ii++) {
              a++; v1 = x[jj*m+ii];
              if (!mxIsNaN(v1)) {c++; v+=v1;}
            }
        if ((a-c)*4<a) *y=v/c; else *y=*x1;
      }
    }
}
