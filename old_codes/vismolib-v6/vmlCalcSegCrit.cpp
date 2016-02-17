
#include <windows.h>
#include <math.h>
#include <stdint.h>
#include "mex.h"

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *sepcrit, *map=NULL;
  uint8_t *data;
  uint8_t d1, d2;
  int *ry;
  double v,sig;
  int nz,r,m,n,x,y,xx,yy,xmin,xmax,ymin,ymax,xxmin,xxmax,yymin,yymax,i,ilo,ihi,pad,j1;
  bool in1;
  
  if (!((nrhs==6 && nlhs<=1) || (nrhs==7 && nlhs==2))) 
    mexErrMsgTxt("call  [hist,sepcrit[,map]] = vmlAnalyzeTile4Seg(data,zzh,xlim,ylim,r,halfgap[,j1])");
  for (i=1;i<nrhs;i++) if (!mxIsDouble(prhs[i])) mexErrMsgTxt("all input parameters must be double");
  for (i=2;i<4;i++) if (mxGetNumberOfElements(prhs[i])!=2) mexErrMsgTxt("need 2 entries in xlim and ylim");
  
  if (mxGetClassID(prhs[0])!=mxUINT8_CLASS) mexErrMsgTxt("data must be of class uint8");
  m = (int)mxGetM(prhs[0]);
  n = (int)mxGetN(prhs[0]);
  data = (uint8_t*)mxGetData(prhs[0]);
  
  nz = (int)mxGetScalar(prhs[1]);
  
  xmin = (int)mxGetPr(prhs[2])[0]-1; xmax = (int)mxGetPr(prhs[2])[1];
  ymin = (int)mxGetPr(prhs[3])[0]-1; ymax = (int)mxGetPr(prhs[3])[1];
  
  v = mxGetScalar(prhs[4]); r = (int)floor(v);
  v = v*v;
  ry = (int*)mxMalloc((r+1)*sizeof(int));
  ry[0]=r;
  for (i=1;i<=r;i++) ry[i] = (int)floor(sqrt(v-i*i));
  
  pad = (int)mxGetScalar(prhs[5]);
  
  plhs[0] = mxCreateDoubleMatrix(1,nz,mxREAL);
  sepcrit = mxGetPr(plhs[0]);

  if (nrhs>=7) {
    j1 = (int)mxGetScalar(prhs[6])-1;
    plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
    map = mxGetPr(plhs[1]);
  }
  
  for(x=0;x<n;x++) 
    for (y=0;y<m;y++)
      if ((d1=data[x*m+y])>0) {
        in1 = (x>=xmin && x<xmax && y>=ymin && y<ymax);
        if (in1) {xxmin=x; xxmax=min(n,x+r+1);} 
        else {xxmin = max(x+(x<xmin || x>=xmax || y>=ymax),xmin); xxmax=min(xmax,x+r+1);}
        for (xx=xxmin;xx<xxmax;xx++) {
          if (in1) {
            if (xx==x) {yymin=y+1; yymax=min(m,y+r+1);}
            else {yymin=max(0,y-ry[xx-x]); yymax=min(m,y+ry[xx-x]+1);}
          } else {
            if (xx==x) {yymin=ymin; yymax=min(ymax,y+r+1);}
            else {yymin=max(ymin,y-ry[xx-x]); yymax=min(ymax,y+ry[xx-x]+1);}
          }
          for (yy=yymin;yy<yymax;yy++) 
            if ((d2=data[xx*m+yy])>0) {
              if (d1<d2) {ilo=(int)d1+pad; ihi=(int)d2-1-pad; sig=1.;} 
              else {ilo=(int)d2+pad; ihi=(int)d1-1-pad; sig=-1.;}
              if (ilo<ihi) {
                if (!in1 || !(xx>=xmin && xx<xmax && yy>=ymin && yy<ymax)) v=0.5; else v=1.;
                for (i=ilo;i<ihi;i++) sepcrit[i]+=v;
                if (map && j1>=ilo && j1<ihi) {
                  map[x*m+y]-=sig*v; map[xx*m+yy]+=sig*v;
                }
              }
            }
        }
      }
      
  mxFree(ry);
}
