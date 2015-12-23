
#include <math.h>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *px, *jorig, *j, *jret;
  double dx1, dx2, x1, x2, d, d1;
  int m, n, mn, nj, i, o, p, p1, p2, pnew, single_step, improved;
  
  if (nrhs!=6) mexErrMsgTxt("usage: j = closest_in_plane(plane_px,jorig,j,dx,dim,single_step)");
  if (nlhs>1) mexErrMsgTxt("usage: j = closest_in_plane(plane_px,jorig,j,dx,dim,single_step)");
  
  if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4])!=2) mexErrMsgTxt("dim must contain two values");
  m = (int)mxGetPr(prhs[4])[0];
  n = (int)mxGetPr(prhs[4])[1];
  mn = m*n;
  
  if (!mxIsDouble(prhs[0]) || mxGetM(prhs[0])!=mn || mxGetN(prhs[0])!=2) mexErrMsgTxt("plane_px must be a prod(dim) x 2 matrix");
  px = mxGetPr(prhs[0]);
  
  if (!mxIsDouble(prhs[1])) mexErrMsgTxt("jorig must be a vector");
  jorig = mxGetPr(prhs[1]);
  nj = mxGetNumberOfElements(prhs[1]);
  
  if (!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2])!=nj) mexErrMsgTxt("j must be a vector of same length as jorig");
  j = mxGetPr(prhs[2]);
  
  if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3])!=2) mexErrMsgTxt("dx must contain two values");
  dx1 = mxGetPr(prhs[3])[0];
  dx2 = mxGetPr(prhs[3])[1];
  
  if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5])!=1) mexErrMsgTxt("single_step must be a scalar");
  single_step = (int)mxGetScalar(prhs[5]);
  
  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[1]),mxGetN(prhs[1]),mxREAL);
  jret = mxGetPr(plhs[0]);
  
  for (i=0;i<nj;i++) {
    o = (int)(*(jorig++))-1;
    if (o<0 | o>=mn) mexErrMsgTxt("invalid value found in jorig");
    x1 = px[o]+dx1; x2 = px[o+mn]+dx2;
    p = (int)(*(j++))-1;
    if (p<0 | p>=mn) mexErrMsgTxt("invalid value found in j");
    d = (px[p]-x1)*(px[p]-x1)+(px[p+mn]-x2)*(px[p+mn]-x2);
    
    for (improved=1;improved;) {
      p1 = p % m; p2 = p / m;
      improved = 0;
      if (p1>0) {
        pnew = p-1;
        d1 = (px[pnew]-x1)*(px[pnew]-x1)+(px[pnew+mn]-x2)*(px[pnew+mn]-x2);
        if (d1<d) {p=pnew; d=d1; improved=1;}
      }
      if (p1<m-1) {
        pnew = p+1;
        d1 = (px[pnew]-x1)*(px[pnew]-x1)+(px[pnew+mn]-x2)*(px[pnew+mn]-x2);
        if (d1<d) {p=pnew; d=d1; improved=1;}
      }
      if (p2>0) {
        pnew = p-m;
        d1 = (px[pnew]-x1)*(px[pnew]-x1)+(px[pnew+mn]-x2)*(px[pnew+mn]-x2);
        if (d1<d) {p=pnew; d=d1; improved=1;}
      }
      if (p2<n-1) {
        pnew = p+m;
        d1 = (px[pnew]-x1)*(px[pnew]-x1)+(px[pnew+mn]-x2)*(px[pnew+mn]-x2);
        if (d1<d) {p=pnew; d=d1; improved=1;}
      }
      if (single_step) break;
    }
    *(jret++) = p+1;
  }
}
