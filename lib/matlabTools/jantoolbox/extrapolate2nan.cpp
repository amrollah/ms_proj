
#include <windows.h>
#include <math.h>
#include "mex.h"

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x, *y, *x1;
  double v, v1, w, sqrt05;
  int m,n,i,j,ii,jj;
  
  if (nrhs!=1 || nlhs>1 || !mxIsDouble(prhs[0])) mexErrMsgTxt("call  y = extrapolate2nan(x)");
  
  x = mxGetPr(prhs[0]);
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
  y = mxGetPr(plhs[0]);
  
  sqrt05 = sqrt(0.5);
  for(x1=x,j=0;j<n;j++) 
    for (i=0;i<m;i++,x1++,y++) {
      if (!mxIsNaN(*x1)) *y=*x1; else {
        for (w=0.,v=0.,jj=max(0,j-1);jj<min(n,j+2);jj++)
          for (ii=max(0,i-1);ii<min(m,i+2);ii++) {
              v1 = x[jj*m+ii];
              if (!mxIsNaN(v1))
                if (ii==i || jj==j) {w+=1.; v+=v1;} 
                else {w+=sqrt05; v+=sqrt05*v1;} 
            }
        if (w>0.) *y=v/w; else *y=*x1;
      }
    }
}
