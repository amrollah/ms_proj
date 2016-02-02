/*
function b = enforcemindistance(x,d,drel)
% b = enforcemindistance(x,d,drel)
%  Returns a logical vector b such that b(1) and b(end) are true and b(i) is true
%  iff x(i)>=x(j)+d and x(i)>=x(j)*(1+drel) where j is the previous true index.
*/

#include <math.h>
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[] )
{
  int i, j=0, nx;
  double *x;
  double d, drel, prev;
  mxLogical *b, *bb;
  
  if (nrhs!=3 || nlhs>1) mexErrMsgTxt("use b = enforcemindistance(x,d,drel)");
    
  x = mxGetPr(prhs[0]);  
  nx = mxGetNumberOfElements(prhs[0]);
  d = mxGetScalar(prhs[1]);  
  drel = 1.0+mxGetScalar(prhs[2]);
      
  plhs[0] = mxCreateLogicalMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]));
  bb = mxGetLogicals(plhs[0]);
  b = bb;
  *(b++) = 1; //first element
  prev = *(x++);
  
  for (i=1; i<nx-1; i++, x++, b++)
  {
    if (*x>=prev+d && *x>=prev*drel) {prev = *x; *b = 1; j = i;}
  }
  *b = 1; //last element
  if ((*x<prev+d || *x<prev*drel) && j>0) bb[j]=0;
}
