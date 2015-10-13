
#include <math.h>
#include <windows.h>
#include "mex.h"

using namespace std;

#define MAX_OUT_FRAC 0.2
#define MIN_CONF_OUT 1.28
#define TARGET_REDUCT_OUT 0.5
#define TRACK_TARGET_REDUCT_OUT 0.9
#define MIN_CONF_GAP 1.28
#define TARGET_REDUCT_GAP 0.1
#define TRACK_TARGET_REDUCT_GAP 0.5

double *hist;
double target_reduct;

typedef struct { 
  int l, r;
  double m, m2, s; 
} range;

void init_zero_range(range &g, int pos) {
  g.l = pos; g.r = pos;
  g.m = 0.; g.m2 = 0.; g.s = 0.;
}

void extend_range_to_right(range &g) {
  double x = hist[g.r];
  g.m += x; g.m2 += x*x; g.s += x*g.r;
  g.r++;
}

void extend_range_to_left(range &g) {
  g.l--;
  double x = hist[g.l];
  g.m += x; g.m2 += x*x; g.s += x*g.l;
}

void shrink_range_from_left(range &g) {  
  double x = hist[g.l];
  g.m -= x; 
  if (g.m<=0.) {g.m2 = 0.; g.s = 0.;}
  else {
    g.m2 -= x*x; if (g.m2<0.) g.m2=0.;
    g.s -= x*g.l; if (g.s<0.) g.s=0.;
  }
  g.l++;
}

void compute_right_range(range &both, range &l, range &r) {
  r.r = both.r; r.l = l.r;
  r.m = both.m-l.m; r.m2 = both.m2-l.m2; r.s = both.s-l.s;
}

void compute_left_range(range &both, range &l, range &r) {
  l.l = both.l; l.r = r.l;
  l.m = both.m-r.m; l.m2 = both.m2-r.m2; l.s = both.s-r.s;
}

void sum_range(range &sum, range &l, range &r) {
  if (r.l!=l.r) mexErrMsgTxt("sum range: left and right part are not adjacent unexpectedly");
  sum.l = l.l; sum.r = r.r;
  sum.m = l.m+r.m; sum.m2 = l.m2+r.m2; sum.s = l.s+r.s;
}

void copy_range(range &g, range &gsrc) {
  g.l = gsrc.l; g.r = gsrc.r;
  g.m = gsrc.m; g.m2 = gsrc.m2; g.s = gsrc.s;
}

double fscore_out(range &out, range &main) {
  //assign and check widths
  int wout = out.r-out.l, wmain = main.r-main.l;
  if (wout<2 || wmain<2) return -2000.0;
  
  //average masses of the ranges
  double avgout = out.m/wout, avgmain = main.m/wmain;
  if (avgout>=avgmain) return -1000.0;
  
  //confidence margins based on normal distribution statistics
  double margin = out.m2-wout*avgout*avgout;
  if (margin<1.0) margin = 1.0; else margin = sqrt(margin);
  double margin_tmp = main.m2-wmain*avgmain*avgmain;
  if (margin_tmp<1.0) margin_tmp = 1.0; else margin_tmp = sqrt(margin_tmp);
  margin = margin/(wout-1)+margin_tmp/(wmain-1);
  
  //prio 1: confidence level high enough
  double conf = (avgmain-avgout)/margin;
  if (conf<MIN_CONF_OUT) return conf-MIN_CONF_OUT;
  
  //prio 2: confident reduction in average mass
  double reduct = 1.0-(avgout+MIN_CONF_OUT*margin)/avgmain;
  if (reduct<target_reduct) return reduct;
  double target_reduct_new = TARGET_REDUCT_OUT+(reduct-TARGET_REDUCT_OUT)*TRACK_TARGET_REDUCT_OUT;
  if (target_reduct_new>target_reduct) target_reduct=target_reduct_new;
  
  //prio 3: maximize width of out range
  return reduct+wout;
}

double fscore_gap(range &gap, range &peak1, range &peak2) {
  //assign and check widths
  int wgap = gap.r-gap.l, w1 = peak1.r-peak1.l, w2 = peak2.r-peak2.l;
  if (wgap<2 || w1<2 || w2<2) return -2000.0;
  
  //average masses of the ranges
  double avggap = gap.m/wgap, avg1 = peak1.m/w1, avg2 = peak2.m/w2;
  if (avggap>=avg1 || avggap>=avg2) return -1000.0;
  
  //confidence margins based on normal distribution statistics
  double margin1 = peak1.m2-w1*avg1*avg1;
  if (margin1<1.0) margin1 = 1.0; else margin1 = sqrt(margin1);
  double margin_gap = gap.m2-wgap*avggap*avggap;
  if (margin_gap<1.0) margin_gap = 1.0; else margin_gap = sqrt(margin_gap);
  margin_gap /= (wgap-1);
  margin1 = margin1/(w1-1)+margin_gap;
  double margin2 = peak2.m2-w2*avg2*avg2;
  if (margin2<1.0) margin2 = 1.0; else margin2 = sqrt(margin2);
  margin2 = margin2/(w2-1)+margin_gap;
  
  //prio 1: confidence level high enough
  double conf = (avg1-avggap)/margin1, conf2 = (avg2-avggap)/margin2;
  if (conf2<conf) conf=conf2;
  if (conf<MIN_CONF_GAP) return conf-MIN_CONF_GAP;
  
  //prio 2: confident reduction in average mass
  double reduct = 1.0-(avggap+MIN_CONF_OUT*margin1)/avg1;
  double reduct2 = 1.0-(avggap+MIN_CONF_OUT*margin2)/avg2;
  if (reduct2<reduct) reduct=reduct2;
  if (reduct<target_reduct) return reduct;
  double target_reduct_new = TARGET_REDUCT_OUT+(reduct-TARGET_REDUCT_OUT)*TRACK_TARGET_REDUCT_GAP;
  if (target_reduct_new>target_reduct) target_reduct=target_reduct_new;
  
  //prio 3: maximize width of gap range
  return reduct+wgap;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *y, *res;
  double y_resolution, ymin, ymax, p, q, x, score, bestscore;
  range all, left, right, left_out, right_out, left_main, right_main;
  range left_peak, right_peak, ltry, rtry, gapltry, gaprtry, gaptry;
  int n, nbins, bin0, i, j, gap_l, gap_r, w;
  bool single_peak;
  
  if (nrhs!=3 || nlhs>2) mexErrMsgTxt("call [res, hist] = otsu_gap(yy, y_resolution, wmin)");
  
  if (!mxIsDouble(prhs[0])) mexErrMsgTxt("yy must contain doubles");
  if (!mxIsDouble(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1) mexErrMsgTxt("y_resolution be scalar");
  if (!mxIsDouble(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1) mexErrMsgTxt("wmin be scalar");
  
  y = mxGetPr(prhs[0]);
  n = mxGetNumberOfElements(prhs[0]);
  y_resolution = mxGetScalar(prhs[1]);
  w = (int)mxGetScalar(prhs[2]);
  
  plhs[0] = mxCreateDoubleMatrix(7,2,mxREAL);
  res = mxGetPr(plhs[0]);
    
  //find data minimum and maximum
  ymin = y[0]; ymax = y[0];
  for (i=1;i<n;i++) {
    if (y[i]<ymin) ymin=y[i];
    if (y[i]>ymax) ymax=y[i];
  }
  res[5] = ymin; res[12] = ymax; //absolute minimum and maximum
  ymin = floor(.5+ymin/y_resolution);
  ymax = floor(.5+ymax/y_resolution);
  bin0 = (int)ymin;
  nbins = 1+(int)ymax-bin0;
  
  //compute data histogram
  if (nlhs<2) hist = (double*)mxMalloc(nbins*sizeof(double));
  else {
    plhs[1] = mxCreateDoubleMatrix(nbins,2,mxREAL);
    hist = mxGetPr(plhs[1]);
  }
  for (i=0;i<nbins;i++) hist[i]=0;
  for (i=0;i<n;i++) {
    p = y[i]/y_resolution;
    q = floor(.5+p);
    j = (int)q-bin0;
    if (p>q) {
      if (j==nbins-1) hist[j]+=1.;
      else {hist[j+1]+=p-q; hist[j]+=1.-p+q;}
    } else {
      if (j==0) hist[j]+=1.;
      else {hist[j-1]+=q-p; hist[j]+=1.-q+p;}
    }
  }
  
  //initialize "all" range info
  all.l = 0; all.r = nbins; all.m = n; all.m2 = 0.; all.s = 0.;
  for (i=0; i<nbins; ++i) { x = hist[i]; all.s += x*i; all.m2 += x*x; }
  
  //if nbins<wmin+4, go right away for a single peak
  single_peak = (nbins<w+4);
  
  //Otsu's method for separating two peaks
  if (!single_peak) {
    init_zero_range(left, 0);
    bestscore = 0.; init_zero_range(ltry, 0);
    for (ltry.r=0; ltry.r<all.r-2;) {
      extend_range_to_right(ltry);
      if (ltry.m<0.5 || ltry.r<2) continue;
      if (ltry.m+0.5>all.m) break;
      compute_right_range(all,ltry,rtry);
      score = ltry.s/ltry.m - rtry.s/rtry.m;
      score = ltry.m * rtry.m * score * score;
      if (score>bestscore) {bestscore=score; copy_range(left,ltry);}
    }
    compute_right_range(all,left,right);
    if (left.r==0) single_peak = true;
  }
  
  init_zero_range(left_out, 0);
  init_zero_range(right_out, nbins);
  
  if (!single_peak) {
    //clean left low density area of left part
    if (left.r-left.l>=4) {
      bestscore = 0.; init_zero_range(ltry, 0);
      target_reduct = TARGET_REDUCT_OUT;
      for (ltry.r=0; ltry.r<left.r-2;) {
        extend_range_to_right(ltry);
        if (ltry.m>MAX_OUT_FRAC*all.m || ltry.m>=0.5*left.m) break;
        compute_right_range(left,ltry,rtry);
        score = fscore_out(ltry,rtry);
        if (score>bestscore) {bestscore=score; copy_range(left_out,ltry);}
      }
    }
    compute_right_range(left,left_out,left_main);
    
    //clean right low density area of right part
    if (right.r-right.l>=4) {
      bestscore = 0.; init_zero_range(rtry, nbins);
      target_reduct = TARGET_REDUCT_OUT;
      for (rtry.l=nbins; rtry.l>right.l+2;) {
        extend_range_to_left(rtry);
        if (rtry.m>MAX_OUT_FRAC*all.m|| rtry.m>=0.5*right.m) break;
        compute_left_range(right,ltry,rtry);
        score = fscore_out(rtry,ltry);
        if (score>bestscore) {bestscore=score; copy_range(right_out,rtry);}
      }
    }
    compute_left_range(right,right_main,right_out);
    
    //left peak comprises at least left 50% of left main, but then until a local maximum of average mass
    bestscore = 0.; gap_l = 0; init_zero_range(ltry, left_main.l);
    for (ltry.r=left_main.l; ltry.r<left_main.r;) {
      extend_range_to_right(ltry);
      if (ltry.m<0.5*left_main.m || ltry.r<left_main.l+2) continue;
      score = ltry.m/(ltry.r-ltry.l);
      if (score>bestscore) {bestscore=score; gap_l=ltry.r;}
    }
    if (gap_l==0) mexErrMsgTxt("computing the left peak failed unexpectedly");
    
    //right peak comprises at least right 50% of right main, but then until a local maximum of average mass
    bestscore = 0.; gap_r = nbins; init_zero_range(rtry, right_main.r);
    for (rtry.l=right_main.r; rtry.l>right_main.l;) {
      extend_range_to_left(rtry);
      if (rtry.m<0.5*right_main.m || rtry.l>right_main.r-2) continue;
      score = rtry.m/(rtry.r-rtry.l);
      if (score>bestscore) {bestscore=score; gap_r=rtry.l;}
    }
    if (gap_r==nbins) mexErrMsgTxt("computing the right peak failed unexpectedly");
    
    //finally optimize the gap
    copy_range(left_peak,left_main); copy_range(right_peak,right_main);
    bestscore = 0.;
    target_reduct = TARGET_REDUCT_GAP;
    for (; w<=gap_r-gap_l; w++) {
      //find best score with required width
      init_zero_range(gapltry,max(left_main.r-w,gap_l));
      init_zero_range(gaprtry,left_main.r);
      for (i=0;i<w;i++) {
        if (gapltry.r<left_main.r) extend_range_to_right(gapltry);
        else extend_range_to_right(gaprtry);
      }
      while (true) {
        compute_left_range(left_main,ltry,gapltry);
        compute_right_range(right_main,gaprtry,rtry);
        sum_range(gaptry,gapltry,gaprtry);
        score = fscore_gap(gaptry,ltry,rtry);
        if (score>bestscore) {
          bestscore=score; 
          copy_range(left_peak,ltry); copy_range(right_peak,rtry);
          x = gaptry.m/w; //average gap mass
        }
        if (gapltry.l>=gapltry.r || gaprtry.r>=gap_r) break;
        shrink_range_from_left(gapltry);
        extend_range_to_right(gaprtry);
      }
    }
    
    //write back results
    p = y_resolution*(left_peak.s/left_peak.m+bin0); //center of left peak
    q = y_resolution*(right_peak.s/right_peak.m+bin0); //center of right peak
    
    res[6] = y_resolution*(-0.5+left.r+bin0); //result of Otsu
    res[4] = y_resolution*(-0.5+left_out.r+bin0); //cleaned minimum
    res[11] = y_resolution*(-0.5+right_out.l+bin0); //cleaned maximum
    
    if (left_peak.r<right_peak.l) { //gap found
      res[0] = y_resolution*(-0.5+right_peak.l+bin0);
      res[7] = y_resolution*(-0.5+left_peak.r+bin0);
      res[2] = q; res[9] = p;
      res[3] = x; res[10] = x; //average gap mass
      res[1] = (right_peak.m/(right_peak.r-right_peak.l)-x)/(q-res[0]);
      res[8] = (left_peak.m/(left_peak.r-left_peak.l)-x)/(res[7]-p);
    } else { //no gap
      res[0]=res[4]; res[7]=res[11];
      res[2] = p; res[9] = q;
      if (left_out.m==0.) res[3] = 0.; else res[3] = left_out.m/(left_out.r-left_out.l);
      if (right_out.m==0.) res[10] = 0.; else res[10] = right_out.m/(right_out.r-right_out.l);
      res[1] = (left_peak.m/(left_peak.r-left_peak.l)-res[3])/(p-res[0]);
      res[8] = (right_peak.m/(right_peak.r-right_peak.l)-res[10])/(res[7]-q);
    }
  } else { //single peak from the beginning
    p = y_resolution*(all.s/all.m+bin0); //center of data
    q = all.m/(all.r-all.l); //average mass
    res[6] = p;
    res[4] = y_resolution*(-0.5+bin0);
    res[11] = y_resolution*(-0.5+nbins+bin0);
    res[0]=res[4]; res[7]=res[11];
    res[2] = p; res[9] = p;
    res[3] = 0.; res[10] = 0.;
    res[1] = q/(p-res[0]);
    res[8] = q/(res[7]-p);
  }
    
  //reshape data histogram to output format if desired
  if (nlhs<2) mxFree(hist);
  else {
    for (i=0;i<nbins;i++) {
      hist[nbins+i] = hist[i];
      hist[i] = y_resolution*(i+bin0);
    }
  }
}
