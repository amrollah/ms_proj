
#include <windows.h>
#include <math.h>
#include "mex.h"

#include <stack>

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  mxLogical *x, *y;
  int i,j,p,m,n;
  stack <int> buf;
  
  if (nrhs!=3 || nlhs>1) mexErrMsgTxt("call  y = floodfill(x,istart,jstart)");
  if (!mxIsLogical(prhs[0])) mexErrMsgTxt("x must be logical");
  for (i=1;i<nrhs;i++) if (!mxIsDouble(prhs[i]) || mxGetM(prhs[i])!=1 || mxGetN(prhs[i])!=1) mexErrMsgTxt("istart and jstart must be scalar");
  
  x = mxGetLogicals(prhs[0]);
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  
  i = (int)mxGetScalar(prhs[1])-1;
  j = (int)mxGetScalar(prhs[2])-1;
  
  if (i<0 || i>=m) mexErrMsgTxt("index istart out of range");
  if (j<0 || j>=n) mexErrMsgTxt("index jstart out of range");
  
  plhs[0] = mxCreateLogicalMatrix(m,n);
  y = mxGetLogicals(plhs[0]);
  
  buf.push(i); buf.push(j);

  //simple 4-way flood fill
  while (!buf.empty()) {
    i = buf.top(); buf.pop();
    j = buf.top(); buf.pop();
    p = j*m+i;
    if (x[p] && !y[p]) {
      y[p] = true;
      if (i>0) {buf.push(i-1); buf.push(j);}
      if (i<m-1) {buf.push(i+1); buf.push(j);}
      if (j>0) {buf.push(i); buf.push(j-1);}
      if (j<n-1) {buf.push(i); buf.push(j+1);}
    }
  }
}
