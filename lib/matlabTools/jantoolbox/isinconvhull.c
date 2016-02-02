/*
function [in,dist] = isinconvhull(xt,x,k,nv)
%test if points xt are in the convex hull of x
%with facets k and outer normal vectors nv

d0 = sum(nv.*x(:,k(:,1)),1);
ri = zeros(1,size(k,1));
for i=1:size(xt,2)
  d = max(sum(nv.*xt(:,i+ri),1)-d0);
  if d<1e-7, in(i) = 1; else in(i) = 0; end;
  dist(i) = d;
end;
in = logical(in);
*/

#include <math.h>
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[] )
{
  int i, j, k, K, n, nf, d, dd, d2, dsq;
  double *x, *xt, *facets, *nv, *d0, *px, *pnv, *in, *D;
  double u,v;
  
  if (nrhs>4) mexErrMsgTxt("Too many input arguments");
  if (nlhs>2) mexErrMsgTxt("Too many output arguments");
  
  if(nrhs<1 || !mxIsDouble(prhs[0]) || mxIsEmpty(prhs[0]))
    mexErrMsgTxt("Parameter xt must contain points");
  
  xt = mxGetPr(prhs[0]);  
  n = mxGetN(prhs[0]);
  d = mxGetM(prhs[0]);
  
  if(nrhs<2 || !mxIsDouble(prhs[1]) || mxGetM(prhs[1])!=d)
    mexErrMsgTxt("Parameter x must have correct dimension");
    
  x = mxGetPr(prhs[1]);
  K = mxGetN(prhs[1]);
  
  if(nrhs<3 || !mxIsDouble(prhs[2]) || mxGetN(prhs[2])!=d)
    mexErrMsgTxt("Parameter k must have correct size");
    
  facets = mxGetPr(prhs[2]);
  nf = mxGetM(prhs[2]);
  
  if(nrhs<4 || !mxIsDouble(prhs[3]) || mxGetM(prhs[3])!=d || mxGetN(prhs[3])!=nf)
    mexErrMsgTxt("Parameter nv must have correct size");
  
  nv = mxGetPr(prhs[3]);
  
  d0 = (double*)mxMalloc(nf*sizeof(double));    
  pnv = nv;
  for (i=0;i<nf;i++) {
    k = ((int)facets[i])-1;
    if (k>=K) mexErrMsgTxt("To large index in k");
    px = x+k*d;
    u = 0.0;
    for (j=0;j<d;j++) u+= (*(px++))*(*(pnv++));
    d0[i] = u;
  }
  
  plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
  //mxSetLogical(plhs[0]);
  in = mxGetPr(plhs[0]);
    
  if (nlhs>1) {
    plhs[1] = mxCreateDoubleMatrix(1,n,mxREAL);
    D = mxGetPr(plhs[1]);
  }  
  
  px = xt;
  for (k=0;k<n;k++) {
    pnv = nv;
    v = 0.0;
    for (i=0;i<nf;i++) {
      u = -d0[i];
      for (j=0;j<d;j++) u+= px[j]*(*(pnv++));
      if (u>=1e-7) {
        if (nlhs<=1) {v=1.0; break;}
        if (u>v) v=u;
      }
    }
    if (v==0.0) (*(in++))=1; else in++;
    if (nlhs>1) if (v==0.0) D++; else (*(D++))=v;
    px += d;
  }
  
  mxFree(d0);
}
