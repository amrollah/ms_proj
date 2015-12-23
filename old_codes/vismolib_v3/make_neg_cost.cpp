
#include <math.h>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *y, *res;
  double offset, score;
  int sum, suml, sum1, s1, s2;
  int n, nbins, i, j, improved, p=0;
  int *hist, *ii, *jj, *ssum, *ssuml;
  
  if (nrhs!=3 || nlhs>1) mexErrMsgTxt("call res = make_neg_cost(y,nbins,offset)");
  
  if (!mxIsDouble(prhs[0])) mexErrMsgTxt("y must contain doubles");
  if (!mxIsDouble(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1) mexErrMsgTxt("nbins must be a scalar");
  if (!mxIsDouble(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1) mexErrMsgTxt("offset must be a scalar");
  
  y = mxGetPr(prhs[0]);
  n = mxGetNumberOfElements(prhs[0]);
  nbins = (int)mxGetScalar(prhs[1]);
  offset = mxGetScalar(prhs[2]);
  
  plhs[0] = mxCreateDoubleMatrix(1,2,mxREAL);
  res = mxGetPr(plhs[0]);
  
  hist = new int[nbins];
  ii = new int[nbins];
  jj = new int[nbins];
  ssum = new int[nbins];
  ssuml = new int[nbins];
  
  for (i=0;i<nbins;i++) hist[i]=0;
  for (i=0;i<n;i++) {
    j = (int)(y[i]*nbins);
    if (j>=nbins) j=nbins-1;
    hist[j]++;
  }
  for (i=0;i<nbins;i++) mexPrintf("%i\n",hist[i]);
  
  i=0; j=1; sum = hist[0]; suml = 0;
  score = (double)sum+offset;
  while(1) {
    improved = 0;
    while (j<nbins) {
      sum1=sum+hist[j]; 
      s1=((double)sum1+offset)/((double)(j+1-i))
      if (s1>score) break;
      score = s1;
      j++;
      improved=1;
      sum = sum1;
    }
    if (!improved) {
      while (i<j-1) {
        sum1=sum-hist[i]; 
        s1=((double)sum1+offset)/((double)(j-i-1))
        if (s1>=score) break;
        score = s1;
        suml += hist[i++];
        improved=1;
        sum = sum1;
      }
    }
    if (!improved) {
      ii[p] = i;
      jj[p] = j;
      ssum[p] = sum;
      ssuml[p] = suml;
      p++;
      if ((i=j+1)>=nbins) break;
      j = i+1; sum = hist[i]; suml = 0;
      score = sum+offset;
    }
  }
  
  j = p-1; i = 0;
  sum = n-suml[0];
  score=((double)sum+offset)/((double)(jj[j]-ii[i]))
  while (j>i) {
    
  }
  
  delete[] hist;
  delete[] ii;
  delete[] jj;
  delete[] ssum;
  delete[] ssuml;
}
