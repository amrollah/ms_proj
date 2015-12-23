VMLCONF = [];
VMLCONF.pad_minutes = 10;
VMLCONF.show_dt = 30;
VMLCONF.img_ext = '.jpeg';
VMLCONF.basefolder = 'D:\l2\data\Cavriglia\img\';
VMLCONF.datafolder = 'X:\proj\COSMOS\CloudTracking\Cavriglia\data\';
VMLCONF.img_offset = [270 250];
%VMLCONF.sky_mask = 'image_mask_2015_01_21.mat';
VMLCONF.sky_mask = 'cavriglia_skymask.png';
VMLCONF.calibration = {'Cavriglia_InternalCalib_20150417.mat',...
  'Cavriglia_Ecalib_20150717.mat',...
  'Cavriglia_model3D.mat'};
VMLCONF.clear_sky_power_model = 'Pclearsky.mat';
VMLCONF.rsun_px = 150;
VMLCONF.gaussblur4sundect_px = 77;
VMLCONF.minangle = 5*pi/180;
VMLCONF.nfilter_mvec = 10;
VMLCONF.nfilter_short_mvec = 5;
VMLCONF.evfilter_short_mvec = 4;
VMLCONF.filter_short_maxstd = 0.2;
VMLCONF.r2b.min = -1;
VMLCONF.r2b.max = 1;
VMLCONF.oflow.mask_width = 2;   %width for optical flow mask [m on plane of height 1]
VMLCONF.oflow.mask_cover_center = 1;   %optical flow mask needs to cover at least center
VMLCONF.oflow.masked = 1;       %for optical flow, use image part where clouds come from
VMLCONF.oflow.ignore_diff = 12; %threshold below which red channel differences are ignored
VMLCONF.oflow.w_upwind = 3;     %additional weight for upwind area
VMLCONF.oflow.median_ksize = 3; %KSize for median filter
VMLCONF.thres.nfilter = 5;
VMLCONF.thres.rneigh = 2;
VMLCONF.thres.nsec = 8;
VMLCONF.thres.rinner = 400;
VMLCONF.thres.router = 600;
VMLCONF.thres.r2b_cutoff = .3;
VMLCONF.thres.r2b_low = 0.02;
VMLCONF.thres.lum_cutoff = 0;
VMLCONF.thres.lum_low = 0.01;
VMLCONF.thres.ngrid = 10;
VMLCONF.thres.wavg.nxy = 7;
VMLCONF.thres.wavg.nlum = 1;
VMLCONF.thres.wavg.sxy = .15;
VMLCONF.thres.wavg.slum = 1;
VMLCONF.thres.wavg.offset = .001;
% VMLCONF.thres.gp.nxy = 7;
% VMLCONF.thres.gp.nlum = 1;
% VMLCONF.thres.gp.sxy = .2;
% VMLCONF.thres.gp.slum = inf;
% VMLCONF.thres.gp.prior = 1;
% VMLCONF.thres.gp.bound = 0.1;
% VMLCONF.thres.gp.ksum_cutoff = 0.05;
VMLCONF.pred.n_feature_maps = 3;
VMLCONF.pred.width0 = 0.2;
VMLCONF.pred.width_v = 1/4;
VMLCONF.pred.tband = 1/3;
VMLCONF.pred.tband0 = 10;
VMLCONF.pred.ngrid = 5;
VMLCONF.pred.nsubgrid = 5;
VMLCONF.pred.clearsky_thres = 0.95;
VMLCONF.pred.covered_thres = 0.4;
VMLCONF.svm.C = 1;
VMLCONF.svm.gradb_train = 0.001;
VMLCONF.svm.gradb_pred = 0.1;
