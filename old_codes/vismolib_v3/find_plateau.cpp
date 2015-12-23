
#include <math.h>
#include <windows.h>
#include <vector>
#include "mex.h"

using namespace std;

#define NEG_NORMINV 1.28f

int n, nbins;

double fscore(int sum, int w) {
  return (1.+NEG_NORMINV*sqrt((1.-(double)w/nbins)/max(1,sum)))*
         (w+sum)/w/(nbins-w+n-sum)*(nbins-w);
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *y, *res;
  double y0, score, score1, bestscore = 0.5;
  int sum, sum1, bestsum;
  int nneigh, i, j0, w, l, lprev, lbest = 0, wbest;
  
  if (nrhs!=5 || nlhs>1) mexErrMsgTxt("call res = find_plateau(y,y0,nbins,wmin,nneigh)");
  
  if (!mxIsDouble(prhs[0])) mexErrMsgTxt("y must contain doubles");
  for (i=1;i<nrhs;i++) if (!mxIsDouble(prhs[i]) || mxGetM(prhs[i])!=1 || mxGetN(prhs[i])!=1) 
    mexErrMsgTxt("all arguments except y must be scalar");
  
  y = mxGetPr(prhs[0]);
  n = mxGetNumberOfElements(prhs[0]);
  y0 = mxGetScalar(prhs[1]);
  nbins = (int)mxGetScalar(prhs[2]);
  w = (int)mxGetScalar(prhs[3]);
  nneigh = (int)mxGetScalar(prhs[4]);
  
  plhs[0] = mxCreateDoubleMatrix(1,4,mxREAL);
  res = mxGetPr(plhs[0]);
  
  res[0]=0; res[1]=1; res[2]=0; res[3]=0;
  bestsum = n;
  wbest = nbins;
  
  vector<int> hist(nbins);
  
  for (i=0;i<nbins;i++) hist[i]=0;
  for (i=0;i<n;i++) {
    l = (int)(y[i]*nbins);
    if (l>=nbins) l=nbins-1;
    hist[l]++;
  }
  j0 = (int)(y0*nbins);
  if (j0>=nbins) j0=nbins-1;
  
  l = max(0,j0+1-w);
  for (i=l, sum=0; i<l+w; i++) sum += hist[i];
  score = fscore(sum,w);
  if (score<=bestscore) {bestscore=score; lbest=l; wbest=w; bestsum=sum;}
  for (i=l+1, sum1=sum; i<=j0 && i+w<nbins; i++) {
    sum1 += (hist[i+w-1]-hist[i-1]);
    score1 = fscore(sum1,w);
    if (score1<score || (score1==score && fabs(i+0.5*(w-1)-j0)<fabs(l+0.5*(w-1)-j0) )) {
      score=score1; sum=sum1; l=i;
      if (score<=bestscore) {bestscore=score; lbest=l; wbest=w; bestsum=sum;}
    }
  }
  
  for (w=w+1;w<nbins;w++) {
    lprev = l;
    l = max(0,lprev-nneigh);
    if (lprev==0) sum += hist[w-1]; else sum += hist[lprev-1];
    for (i=l;i<lprev-1;i++) sum += (hist[i]-hist[i+w]);
    score = fscore(sum,w);
    if (score<=bestscore) {bestscore=score; lbest=l; wbest=w; bestsum=sum;}
    for (i=l+1, sum1=sum; i<lprev+nneigh && i+w<=nbins; i++) {
      sum1 += (hist[i+w-1]-hist[i-1]);
      score1 = fscore(sum1,w);
      if (score1<score || (score1==score && fabs(i+0.5*(w-1)-j0)<fabs(l+0.5*(w-1)-j0) )) {
        score=score1; sum=sum1; l=i;
        if (score<=bestscore) {bestscore=score; lbest=l; wbest=w; bestsum=sum;}
      }
    }
  }
  score = (double)bestsum/wbest;
  if (lbest>0) {
    res[0]=(0.5+lbest)/nbins;
    for (i=lbest-1;i>=0;i--) { 
      score1 = (hist[i]-score)/(lbest-i)*nbins; 
      if (score1>res[2]) res[2]=score1;
    }
  }
  if (lbest+wbest<nbins) {
    res[1]=(lbest+wbest-0.5)/nbins;
    for (i=lbest+wbest;i<nbins;i++) { 
      score1 = (hist[i]-score)/(i-(lbest+wbest-1))*nbins; 
      if (score1>res[3]) res[3]=score1;
    }
  }
}
