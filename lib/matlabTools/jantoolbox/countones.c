
#include <math.h>
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[] )
{
  int i, len, cnt=0;
  mxLogical *x;
  int *y;
  
  if (nrhs!=1) mexErrMsgTxt("Single input argument required");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments");
    
  if(nrhs<1 || !mxIsLogical(prhs[0]) || mxIsEmpty(prhs[0]))
    mexErrMsgTxt("Parameter x must contain booleans");
    
  x = mxGetLogicals(prhs[0]);
  len = mxGetNumberOfElements(prhs[0]);
    
  plhs[0] = mxCreateNumericMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]), mxINT32_CLASS, mxREAL); 
  y = mxGetPr(plhs[0]);
  
  for (i=0;i<len;i++)
    if (*(x++)==0.0) {*(y++)=0.0; cnt=0;} else *(y++)=(++cnt);
}
