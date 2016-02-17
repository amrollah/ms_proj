
#include <windows.h>
#include <math.h>
#include "mex.h"

using namespace std;

#include <vector>
#include <algorithm>

template<class T> double median(vector<T> x, size_t n) 
{
  //size_t n = x.size();
  sort(x.begin(), x.begin()+n);
  if (n&1) return x[n>>1];
  else return 0.5*(x[(n>>1)-1]+x[n>>1]);
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x, *y, *w, *y1, *x2, *xy, *coeff;
  double c0, c1, c0_, c1_, XWtXW00, XWtXW01, XWtXW11, XWtYW0, XWtYW1, wsum, w2, deti, sae, sae0, med, alpha;
  int nx, ny1, N, n, i, niter;
  
  if (nrhs!=3 || nlhs>1) mexErrMsgTxt("call coeff = vmlPointwiseL1Regression(x,y,alpha)");
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) mexErrMsgTxt("x and y must be double");
  if (!mxIsDouble(prhs[2]) || mxGetM(prhs[2])!=1  || mxGetN(prhs[2])!=1) mexErrMsgTxt("alpha must be a scalar");
  
  x = mxGetPr(prhs[0]);
  nx = mxGetNumberOfElements(prhs[0]);
  
  if (mxGetM(prhs[1])!=nx) mexErrMsgTxt("number of rows of y does not match number of elements of x");
  y = mxGetPr(prhs[1]);
  N = mxGetN(prhs[1]);
  
  alpha = mxGetScalar(prhs[2]);
  
  plhs[0] = mxCreateDoubleMatrix(2,N,mxREAL);
  coeff = mxGetPr(plhs[0]);
  
  w = (double*)mxMalloc(nx*sizeof(double));
  x2 = (double*)mxMalloc(nx*sizeof(double));
  xy = (double*)mxMalloc(nx*sizeof(double));
  vector<double> yy(nx);
  for (i=0;i<nx;i++) x2[i]=x[i]*x[i];
  
  for (y1=y,n=0;n<N;n++,y1+=nx,coeff+=2) {
    for (i=0,ny1=0;i<nx;i++) 
      if (mxIsNaN(y1[i])) w[i]=-1.; else {w[i]=1.; xy[i]=x[i]*y1[i]; yy[ny1++]=y1[i];}
    c0_ = 1.; c1_ = 1.; c0 = 0.; c1 = 0.;
    for(niter=0;(fabs(c0-c0_)>1.e-6 || fabs(c1-c1_)>1.e-6) && niter<200;niter++) {
      XWtXW00=0.; XWtXW01=0.; XWtXW11=0.; XWtYW0=0.; XWtYW1=0.;
      for (i=0;i<nx;i++) if (w[i]>=0.) {
        w2 = w[i]*w[i];
        XWtXW00+=w2; XWtXW01+=w2*x[i]; XWtXW11+=w2*x2[i];
        XWtYW0+=w2*y1[i]; XWtYW1+=w2*xy[i];
      }
      deti = 1./(XWtXW00*XWtXW11-XWtXW01*XWtXW01);
      c0_ = c0; c1_ = c1;
      c0 = deti*(XWtXW11*XWtYW0-XWtXW01*XWtYW1);
      c1 = deti*(XWtXW00*XWtYW1-XWtXW01*XWtYW0);
      for (wsum=0.,i=0;i<nx;i++) if (w[i]>=0.) {
        w[i]=sqrt(1./max(1.e-6,fabs(c0+c1*x[i]-y1[i])));
        wsum += w[i]; //maxerr1 = max(maxerr1,err1);
      }
      for (i=0;i<nx;i++) w[i] = w[i]/wsum;
    }
    med = median(yy,ny1);
    for (sae=0.,sae0=0.,i=0;i<nx;i++) if (w[i]>=0.) {
      sae += fabs(c0+c1*x[i]-y1[i]);
      sae0 += fabs(med-y1[i]);
    }
    //mexPrintf("c0=%f c1=%f sae0=%f sae1=%f c1pen=%f\n",c0,c1,sae0,sae,fabs(c1)*fabs(x[nx-1]-x[0]));
    //mexPrintf("eq at %f\n",(sae0-sae)/(fabs(c1)*fabs(x[nx-1]-x[0])));
    if (sae+alpha*fabs(c1)*fabs(x[nx-1]-x[0])>=sae0) { //prefer fit with 0 slope / constant model
      coeff[0] = med; coeff[1] = 0.;
    } else { //prefer linear model
      coeff[0]=c0; coeff[1]=c1; //coeff[2]=mae/ny1; coeff[3]=(mae-maxerr1)/(ny1-1);
    }
  }
  
  mxFree(w);
  mxFree(x2);
  mxFree(xy);
}
