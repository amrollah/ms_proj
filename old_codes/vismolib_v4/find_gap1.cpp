
#include <math.h>
#include <windows.h>
#include "mex.h"

using namespace std;

#define NEG_NORMINV 1.28f

int n, nbins;

double fscore(int sum, int w) {
//  return (1.+NEG_NORMINV*sqrt((1.-(double)w/nbins)/max(1,sum)))*
//         (w+sum)/w/(nbins-w+n-sum)*(nbins-w);
  return (double)sum/w - (double)(n-sum)/(nbins-w);
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *y, *res, *hist;
  double y_resolution, score, bestscore, ymin, ymax, maxgapfrac, avglo;
  int bin0, leftbin, sum, sum1, bestsum, maxsum;
  int i, w, l, lprev, lbest, wbest, j0, wmin;
  
  if (nrhs<3 || nlhs>4 || nlhs>2) mexErrMsgTxt("call res = find_gap(yy, y_resolution, wmin, maxgapfrac [default=0.25])");
  
  if (!mxIsDouble(prhs[0])) mexErrMsgTxt("yy must contain doubles");
  for (i=1;i<nrhs;i++) if (!mxIsDouble(prhs[i]) || mxGetM(prhs[i])!=1 || mxGetN(prhs[i])!=1) 
    mexErrMsgTxt("arguments y_resolution, wmin, maxgapfrac must be scalars");
  
  y = mxGetPr(prhs[0]);
  n = mxGetNumberOfElements(prhs[0]);
  y_resolution = mxGetScalar(prhs[1]);
  wmin = (int)mxGetScalar(prhs[2]);
  if (nrhs<4) maxgapfrac = 0.25; else maxgapfrac = mxGetScalar(prhs[3]);
  maxsum = (int)ceil(n*maxgapfrac);
    
  //find data minimum and maximum
  ymin = y[0]; ymax = y[0];
  for (i=1;i<n;i++) {
    if (y[i]<ymin) ymin=y[i];
    if (y[i]>ymax) ymax=y[i];
  }
  ymin = floor(.5+ymin/y_resolution);
  ymax = floor(.5+ymax/y_resolution);
  bin0 = (int)ymin;
  nbins = 1+(int)ymax-bin0;
  
  //initialize function results and histogram
  plhs[0] = mxCreateDoubleMatrix(2,4,mxREAL);
  res = mxGetPr(plhs[0]);
  for (i=0;i<8;i++) res[i] = 0.;
  
  if (nlhs<2) hist = (double*)mxMalloc(2*nbins*sizeof(double));
  else {
    plhs[1] = mxCreateDoubleMatrix(2,nbins,mxREAL);
    hist = mxGetPr(plhs[1]);
  }

  //compute data histogram
  for (i=0;i<nbins;i++) hist[i]=0;
  for (i=0;i<n;i++) hist[(int)floor(.5+y[i]/y_resolution)-bin0]++;
  
  //find and remove outliers to the right
  lbest = nbins; sum = 0; bestsum = 0; bestscore = -(double)n/nbins;
  for (i=nbins-1; i>0; i--) {
    sum += hist[i];
    if (sum>=maxsum) break;
    score = fscore(sum,nbins-i);
    if (score<=bestscore) {bestscore=score; lbest=i; bestsum = sum;}
  }
  if (lbest<nbins) avglo = (double)bestsum/(nbins-lbest); else avglo = 0.;
  res[2]=y_resolution*(lbest+bin0);
  for (i=lbest-1;i>=0;i--) { 
    score = (hist[i]-avglo)/(lbest-i)/y_resolution;
    if (score>res[3]) res[3]=score;
  }
  nbins = lbest;
  n -= bestsum;
    
  //find and remove outliers to the left
  wbest = 0; sum = 0; bestsum = 0; bestscore = -(double)n/nbins;
  for (i=1; i<nbins; i++) {
    sum += hist[i-1];
    if (sum>=maxsum) break;
    score = fscore(sum,i);
    if (score<=bestscore) {bestscore=score; wbest=i; bestsum = sum;}
  }
  if (wbest>0) avglo = (double)bestsum/wbest; else avglo = 0.;
  res[0]=y_resolution*(wbest-1+bin0);
  for (i=wbest;i<nbins;i++) { 
    score = (hist[i]-avglo)/(i-(wbest-1))/y_resolution; 
    if (score>res[1]) res[1]=score;
  }
  leftbin = wbest; nbins -= wbest;
  n -= bestsum;
  
  if (w<=nbins) {
    //Otsu's method for separating two peaks
    for (i=1, sum=0; i<nbins; ++i) sum += i * hist[leftbin+i];
    bestscore = 0.; sum1=0; w = 0;
    for (i=0;i<nbins;i++) {
      w += hist[leftbin+i];
      if (w==0) continue;
      if (w==n) break;
      sum1 += i*hist[leftbin+i];
      score = (double)sum1/w - (double)(sum-sum1)/(n-w);
      score = (double)w*(n-w)*score*score;
      if (score>bestscore) {j0 = i; bestscore=score; }
    }
    
    mexPrintf("%f\n",y_resolution*(leftbin+j0+bin0));

    //search for gap  
    lbest = 0; wbest = nbins; bestscore = -(double)n/nbins;
    for (w=wmin; w<=nbins; w++) {
      //find best score with required width
      l = max(0,j0+1-w);
      for (i=l, sum=0; i<l+w; i++) sum += hist[leftbin+i];
      score = fscore(sum,w);
      if (score<=bestscore) {bestscore=score; lbest=l; wbest=w; bestsum=sum;}
      for (i=l+1, sum1=sum; i<=j0 && i+w<nbins; i++) {
        sum1 += (hist[leftbin+i+w-1]-hist[leftbin+i-1]);
        if (sum1<sum) {
          sum=sum1; l=i; score=fscore(sum,w); 
          if (score<bestscore && sum<=maxsum) {bestscore=score; lbest=l; wbest=w; bestsum=sum;}
        }
      }
      mexPrintf("%i %i %f %f / %f %f\n",w,sum,score,bestscore,y_resolution*(leftbin+l+bin0),y_resolution*(leftbin+l+w-1+bin0));
      if (sum>maxsum) break;
    }

//possibly more efficient but not exhaustive search
//    for (w=w+1;w<nbins && sum<=maxsum;w++) {
//      lprev = l;
//      l = max(0,lprev-nneigh);
//      if (lprev==0) sum += hist[leftbin+w-1]; else sum += hist[leftbin+lprev-1];
//      for (i=l;i<lprev-1;i++) sum += (hist[leftbin+i]-hist[leftbin+i+w]);
//      score = fscore(sum,w);
//      if (score<bestscore && sum<=maxsum) {bestscore=score; lbest=l; wbest=w; bestsum=sum;}
//      for (i=l+1, sum1=sum; i<lprev+nneigh && i+w<=nbins; i++) {
//        sum1 += (hist[leftbin+i+w-1]-hist[leftbin+i-1]);
//        score1 = fscore(sum1,w);
//        if (score1<score) {
//          score=score1; sum=sum1; l=i;
//          if (score<bestscore && sum<=maxsum) {bestscore=score; lbest=l; wbest=w; bestsum=sum;}
//        }
//      }
//    }    

    if (lbest>0 && lbest+wbest<nbins) { //found gap
      avglo = (double)bestsum/wbest;
      res[4]=y_resolution*(leftbin+lbest+bin0);
      for (i=lbest-1;i>=0;i--) { 
        score = (hist[leftbin+i]-avglo)/(lbest-i)/y_resolution;
        if (score>res[5]) res[5]=score;
      }
      res[6]=y_resolution*(leftbin+lbest+wbest-1+bin0);
      for (i=lbest+wbest;i<nbins;i++) { 
        score = (hist[leftbin+i]-avglo)/(i-(lbest+wbest-1))/y_resolution; 
        if (score>res[7]) res[7]=score;
      }
    }
  }
  
  if (nlhs<2) mxFree(hist);
  else {
    nbins = 1+(int)ymax-bin0;
    for (i=nbins-1;i>=0;i--) {
      hist[(i<<1)+1] = hist[i];
      hist[i<<1] = y_resolution*(i+bin0);
    }
  }
}
