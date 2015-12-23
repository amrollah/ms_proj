
#include <math.h>
#include "mex.h"

using namespace std;

double *pyx, *xx, *yy;
double ycur, xcur, dist, distnew, distx, disty;
int m, n, mmid, nmid, mn, nx, ny, ix, iy, j, jx, jy, jnew, jxnew, jynew, dx, dy, nxmid, nymid, searchdone, searchdone_new, r;
int *J;

int searchx(int deltay1) {
  int jx1 = jx, jy1 = jy+deltay1, j1 = j+deltay1, delta1=0;
  double dist1, distl, distr;
  int ret = 0;
  if (pyx[j1+mn]>=nan_thres) { //first walk to not NaN area
    ret = 1;
    if (jx1<nmid) delta1=1; else delta1=-1;
    while(1) {
      jx1 = jx1+delta1; j1=j1+m*delta1;
      if (jx1==0 || jx1==n-1 || pyx[j1+mn]<nan_thres) break;
    }
    if (pyx[j1+mn]>=nan_thres) return 1;
  }
  distx = pyx[j1+mn]-xcur; disty = pyx[j1]-ycur;
  dist1 = distx*distx+disty*disty;
  if (jx==0 || delta1>0) distl = dist1+1.0; else {
    distx = pyx[j1-m+mn]-xcur; disty = pyx[j1-m]-ycur;
    distl = distx*distx+disty*disty;
  }
  if (jx==n-1 || delta1<0) distr = dist1+1.0; else {
    distx = pyx[j1+m+mn]-xcur; disty = pyx[j1+m]-ycur;
    distr = distx*distx+disty*disty;
  }
  if (distr<distl && distr<=dist1) delta1=1;
  else if (distl<distr && distl<=dist1) {delta1=-1; distr=distl;}
  else delta1 = 0;
  if (delta1!=0) {
    while (distr<=dist1) {
      jx1 = jx1+delta1; j1=j1+m*delta1; dist1=distr;
      if (jx1==0 || jx1==n-1) break;
      distx = pyx[j1+delta1*m+mn]-xcur; disty = pyx[j1+delta1*m]-ycur;
      distr = distx*distx+disty*disty;
    }
  }
  if (dist1<distnew) {
    distnew = dist1; jxnew = jx1; jynew = jy1; jnew = j1; searchdone_new = 1;
  }
  return ret;
}

int searchy(int deltax1) {
  int jx1 = jx+deltax1, jy1 = jy, j1 = j+deltax1*m, delta1=0;
  double dist1, distl, distr;
  int ret = 0;
  if (pyx[j1+mn]>=nan_thres) { //first walk to not NaN area
    ret = 1;
    if (jy1<mmid) delta1=1; else delta1=-1;
    while(1) {
      jy1 = jy1+delta1; j1=j1+delta1;
      if (jy1==0 || jy1==m-1 || pyx[j1+mn]<nan_thres) break;
    }
    if (pyx[j1+mn]>=nan_thres) return 1;
  }
  distx = pyx[j1+mn]-xcur; disty = pyx[j1]-ycur;
  dist1 = distx*distx+disty*disty;
  if (jy==0 || delta1>0) distl = dist1+1.0; else {
    distx = pyx[j1-1+mn]-xcur; disty = pyx[j1-1]-ycur;
    distl = distx*distx+disty*disty;
  }
  if (jy==m-1 || delta1<0) distr = dist1+1.0; else {
    distx = pyx[j1+1+mn]-xcur; disty = pyx[j1+1]-ycur;
    distr = distx*distx+disty*disty;
  }
  if (distr<distl && distr<=dist1) delta1=1;
  else if (distl<distr && distl<=dist1) {delta1=-1; distr=distl;}
  else delta1 = 0;
  if (delta1!=0) {
    while (distr<=dist1) {
      jy1 = jy1+delta1; j1=j1+delta1; dist1=distr;
      if (jy1==0 || jy1==m-1) break;
      distx = pyx[j1+delta1+mn]-xcur; disty = pyx[j1+delta1]-ycur;
      distr = distx*distx+disty*disty;
    }
  }
  if (dist1<distnew) {
    distnew = dist1; jxnew = jx1; jynew = jy1; jnew = j1; searchdone_new = 2;
  }
  return ret;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  int at_border, at_border1;
  if (nrhs!=4 || nlhs>1) mexErrMsgTxt("usage: j = closest_in_plane(prg_yy,prg_xx,plane_pyx,j0)");
      
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) mexErrMsgTxt("prg_yy and prg_yy must be vectors");
  yy = mxGetPr(prhs[0]);
  ny = mxGetNumberOfElements(prhs[0]);
  xx = mxGetPr(prhs[1]);
  nx = mxGetNumberOfElements(prhs[1]);
  
  if (!mxIsDouble(prhs[2]) || mxGetN(prhs[2])!=2) mexErrMsgTxt("plane_pyx must be a k x 2 matrix");
  pyx = mxGetPr(prhs[2]);
  k = mxGetM(prhs[2]);

  if (!mxIsDouble(prhs[3]) || mxGetM(prhs[3])!=ny  || mxGetN(prhs[3])!=nx) mexErrMsgTxt("j0 must be a ny x nx matrix");
  J0 = mxGetPr(prhs[3]);
  
  plhs[0] = mxCreateNumericMatrix(ny, nx, mxINT32_CLASS, mxREAL);
  J  = (int*) mxGetData(plhs[0]);
  iy = nx*ny; for (ix=0;ix<iy;ix++) J[ix] = 0;
  
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      ycur = yy[iy]; xcur = xx[ix];
      j = (int)J0[ix*ny+iy];
      if (j>0) {J[ix*ny+iy] = j; continue;}
      for (r=1;;r++) {
      }
    }
  }
    
          
          j = J[(ix-dx)*ny+iy-dy];
          if (j==0) {jx = nmid; jy = mmid; j = jx*m+jy;}
          else {
            if (j>0) j = j-1; else j=-j-1;
            jy = j % m; jx = j / m;
          }
          
          distx = pyx[j+mn]-xcur; disty = pyx[j]-ycur;
          dist = distx*distx+disty*disty;
          searchdone = 0;
          
          while(1) {
            jxnew = jx; jynew = jy; jnew = j; distnew = dist;
            
            if (searchdone!=1) searchx(0);
            if (searchdone!=2) searchy(0);
            
            for (r=1;r<=search_radius && jnew==j;r++) {
              at_border1 = 0;
              if (jy-r<0) at_border1 = 1; else if (searchx(-r)) at_border1 = 1;
              if (jy+r>m-1) at_border1 = 1; else if (searchx(r)) at_border1 = 1;
              if (jx-r<0) at_border1 = 1; else if (searchy(-r)) at_border1 = 1;
              if (jx+r>n-1) at_border1 = 1; else if (searchy(r)) at_border1 = 1;
              if (r==1) at_border = at_border1;
              if (!at_border && r>=search_radius_inner) break;
            }
            
            if (jnew==j) break;
            jx = jxnew; jy = jynew; j = jnew; dist = distnew; searchdone = searchdone_new;
          }
          
          if (jx>0 && pyx[j-m+mn]<nan_thres) distx=fabs(pyx[j+mn]-pyx[j-m+mn]); else distx=0.;
          if (jx<n-1 && pyx[j+m+mn]<nan_thres) distnew=fabs(pyx[j+mn]-pyx[j+m+mn]); else distnew=0.;
          if (distnew>distx) distx=distnew;
          if (jy>0 && pyx[j-1]<nan_thres) disty=fabs(pyx[j]-pyx[j-1]); else disty=0.;
          if (jy<m-1 && pyx[j+1]<nan_thres) distnew=fabs(pyx[j]-pyx[j+1]); else distnew=0.;
          if (distnew>disty) disty=distnew;
          if (dist<=(distx*distx+disty*disty)) J[ix*ny+iy] = 1+j; else J[ix*ny+iy] = -1-j;
        }
      }
    }
  }
}
