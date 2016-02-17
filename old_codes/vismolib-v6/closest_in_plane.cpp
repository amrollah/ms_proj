
#include <math.h>
#include "mex.h"

//#include "nanoflann/include/nanoflann.hpp"

#define nan_thres 1.e11
#define search_radius 15
#define search_radius_inner 4

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
  if (nrhs!=4 || nlhs>1) mexErrMsgTxt("usage: j = closest_in_plane(prg_yy,prg_xx,plane_pyx,ppx_dim)");
  
  if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3])!=2) mexErrMsgTxt("dim must contain two values");
  m = (int)mxGetPr(prhs[3])[0];
  n = (int)mxGetPr(prhs[3])[1];
  mn = m*n;
  
  if (!mxIsDouble(prhs[2]) || mxGetM(prhs[2])!=mn || mxGetN(prhs[2])!=2) mexErrMsgTxt("plane_pyx must be a prod(dim) x 2 matrix");
  pyx = mxGetPr(prhs[2]);
  
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) mexErrMsgTxt("prg_yy and prg_yy must be vectors");
  yy = mxGetPr(prhs[0]);
  ny = mxGetNumberOfElements(prhs[0]);
  xx = mxGetPr(prhs[1]);
  nx = mxGetNumberOfElements(prhs[1]);
  
  nxmid = nx>>1; nymid = ny>>1;
  nmid = n>>1; mmid = m>>1;
  
  plhs[0] = mxCreateNumericMatrix(ny, nx, mxINT32_CLASS, mxREAL);
  J  = (int*) mxGetData(plhs[0]);
  iy = nx*ny; for (ix=0;ix<iy;ix++) J[ix] = 0;
  
  for (dx=-1;dx<=1;dx+=2) {
    for (dy=-1;dy<=1;dy+=2) {
      for (ix=nxmid+((dx-1)/2); ix>=0 && ix<nx; ix+=dx) {
        for (iy=nymid+((dy-1)/2); iy>=0 && iy<ny; iy+=dy) {
          ycur = yy[iy]; xcur = xx[ix];
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



//
//// attempt to solve the problem with KD trees, but it turned out to be not precise enough
//
//using namespace nanoflann;
//
//double *pyx;
//size_t n;
//
//struct PointCloud
//{
//  // Must return the number of data points
//  inline size_t kdtree_get_point_count() const { return n; }
//
//  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
//  inline double kdtree_distance(const double *p1, const size_t idx_p2,size_t /*size*/) const
//  {
//    const double d0=p1[0]-pyx[idx_p2];
//    const double d1=p1[1]-pyx[n+idx_p2];
//    return d0*d0+d1*d1;
//  }
//
//  // Returns the dim'th component of the idx'th point in the class:
//  // Since this is inlined and the "dim" argument is typically an immediate value, the
//  //  "if/else's" are actually solved at compile time.
//  inline double kdtree_get_pt(const size_t idx, int dim) const
//  {
//    if (dim==0) return pyx[idx];
//    else return pyx[n+idx];
//  }
//
//  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
//  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
//  template <class BBOX>
//  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
//};
//
//
//void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
//{
//  double *xx, *yy;
//  size_t nx, ny, ix, iy, j_out;
//  double dist_out;
//  int *J;
//  double query_point[2];
//  
//  if (nrhs!=3 || nlhs>1) mexErrMsgTxt("usage: j = closest_in_plane(prg_yy,prg_xx,plane_pyx)");
//  
//  if (!mxIsDouble(prhs[2]) || mxGetN(prhs[2])!=2) mexErrMsgTxt("plane_pyx must be a prod(dim) x 2 matrix");
//  pyx = mxGetPr(prhs[2]);
//  n = mxGetM(prhs[2]);
//  
//  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) mexErrMsgTxt("prg_yy and prg_yy must be vectors");
//  yy = mxGetPr(prhs[0]);
//  ny = mxGetNumberOfElements(prhs[0]);
//  xx = mxGetPr(prhs[1]);
//  nx = mxGetNumberOfElements(prhs[1]);
//  
//  plhs[0] = mxCreateNumericMatrix(ny, nx, mxINT32_CLASS, mxREAL);
//  J  = (int*) mxGetData(plhs[0]);
//  
//  PointCloud cloud;
//
//  // construct a kd-tree index:
//  typedef KDTreeSingleIndexAdaptor<
//    L2_Simple_Adaptor<double, PointCloud > , PointCloud, 3 /* dim */
//    > my_kd_tree_t;
//
//  my_kd_tree_t   index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(1 /* max leaf */) );
//  index.buildIndex();
//
//  for (ix=0; ix<nx; ix++) {
//    for (iy=0; iy<ny; iy++) {
//      query_point[0] = yy[iy]; query_point[1] = xx[ix];
//      index.knnSearch(query_point, 1, &j_out, &dist_out);
//      J[ix*ny+iy] = 1+j_out;
//    }
//  }
//
//}
