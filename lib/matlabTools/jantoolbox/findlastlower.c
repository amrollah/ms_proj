/*
function z = findlastlower(x,y,bSorted)
% z = findlastlower(x,y,bSorted)
%  For each x(i), find the last element of y that is smaller
%  than x(i).
%  y must be sorted.
%  If x is sorted, you can set bSorted=1, in this case
%  the algorithm runs in time O(size(x)+size(y)), otherwise
%  O(size(x)*size(y)).
%  After the call size(z)=size(x), and z(i)=j if y(j)<x(i)
%  and y(j+1)>x(i), while z(i)=j+0.5 if y(j)=x(i)
%  and y(j-1)<x(i).
*/

#include <math.h>
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],
	          int nrhs, const mxArray *prhs[] )
{
  int i, j=0, nx, ny, bSorted=0;
  double *x, *y, *z;
  
  if (nrhs>3) mexErrMsgTxt("Too many input arguments");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments");
    
  if(nrhs<1 || !mxIsDouble(prhs[0])) mexErrMsgTxt("x must contain doubles");    
  if(nrhs<2 || !mxIsDouble(prhs[1])) mexErrMsgTxt("y must contain sorted doubles");
  
  if (nrhs==3)
  {
    if(!mxIsDouble(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1) 
      mexErrMsgTxt("bSorted must contain 0 or 1");        
    bSorted = (int)mxGetScalar(prhs[2]);
  }
    
  x = mxGetPr(prhs[0]);  
  nx = mxGetNumberOfElements(prhs[0]);
  y = mxGetPr(prhs[1]);  
  ny = mxGetNumberOfElements(prhs[1]);
      
  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),mxREAL);
  z = mxGetPr(plhs[0]);
  
  for (i=0; i<nx; i++)
  {
    if (!bSorted) j=0;
    while (j<ny && y[j]<x[i]) j++;
    if (j<ny && y[j]==x[i]) z[i]=(double)j+0.5; else z[i]=(double)j;
  }
}
