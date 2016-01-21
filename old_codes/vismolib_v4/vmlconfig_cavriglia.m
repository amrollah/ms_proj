VMLCONF = [];
VMLCONF.pad_minutes = 10;
VMLCONF.show_dt = 30;
VMLCONF.img_ext = '_Debevec.jpeg';
VMLCONF.basefolder = 'U:\HDRImages_Cavriglia\img\';
%VMLCONF.basefolder = 'D:\enc_nobak\data\Cavriglia\img\';
VMLCONF.tmpfolder = 'C:\data\cav\log_data\';
VMLCONF.datafolder = 'C:\Users\chamsei\Documents\GitHub\ms_proj\data\Cavriglia\';
VMLCONF.img_offset = [270 250];
%VMLCONF.sky_mask = 'image_mask_2015_01_21.mat';
VMLCONF.sky_mask = 'cavriglia_skymask.png';
VMLCONF.calibration = {'Cavriglia_InternalCalib_20150417.mat',...
  'Cavriglia_Ecalib_20150717.mat',...
  'Cavriglia_model3D.mat'};
VMLCONF.TimeZone = 'Europe/Zurich';
VMLCONF.plant_coords = [50,-135,0;195,-135,0;340,-135,0;342,-58,0;345,20,0;188,25,0;32,30,0;40,-53,0]';

VMLCONF.rsun_px = 150;
VMLCONF.rsun_px_close = 50;
VMLCONF.minangle = 5*pi/180;
VMLCONF.daylight_thres = 0.05;

VMLCONF.sundetect.blur_px = 77;
VMLCONF.sundetect.blur2_px = 17;
VMLCONF.sundetect.r_roi_px = 50;
VMLCONF.sundetect.r1_px = 15;
VMLCONF.sundetect.r2_px = 20;
VMLCONF.sundetect.star_r1_px = 10;
VMLCONF.sundetect.star_r2_px = 100;
VMLCONF.sundetect.satsun_thres = 50;
VMLCONF.sundetect.star_thres_cov = 0.05;
VMLCONF.sundetect.star_thres_clear = 0.08;
VMLCONF.sundetect.maxerr_px = 12;
VMLCONF.glare_rect_w = 8;
VMLCONF.sundetect.r_mfi_px = 20;

VMLCONF.mot.nframes = 5; %number of recent frames to use
VMLCONF.mot.w_cc_newest = 1/3; %weight of the most recent cloud coverage
VMLCONF.mot.zero_tol = 0;

VMLCONF.oflow.vmax = .10; %max speed [m/sec on the 1m plane]
VMLCONF.oflow.vstep = 0.005; %initial step [m/sec on the 1m plane]
VMLCONF.oflow.xtol_rel = 0.01; %relative tolerance for the optimization
VMLCONF.oflow.heigh_weight_width = 2;   %width of optical flow high weight area [m on plane of height 1]
VMLCONF.oflow.heigh_weight_cover_center = 1;   %optical flow high weight area needs to cover at least center
VMLCONF.oflow.low_weight = 0.2; %low weight value (heigh weight value is 1)

VMLCONF.mfi.nfilter = 5;

VMLCONF.seg.ngrid = 25;
VMLCONF.seg.tile_max_nan = 0.5;
VMLCONF.seg.alpha = 0.1;
VMLCONF.seg.outside_mass = 0.1;
VMLCONF.seg.extreme_outside_mass = 0.02;
VMLCONF.seg.n_relax_outside = 3;
VMLCONF.seg.r2b_midrange = [0 0.5];
VMLCONF.seg.r2b_resolution = 0.3;
VMLCONF.seg.r2b_hist_resolution = 0.05;
VMLCONF.seg.wsmooth = .25;
VMLCONF.seg.w_local_min = .05;
VMLCONF.seg.w_boundary = 8;
VMLCONF.seg.r2b_dmax = 0.2;
VMLCONF.seg.rsun_mfi4outer = 80;
VMLCONF.seg.r2b_multistart_step = .3;
VMLCONF.seg.max_sunglare_remove = 0.1; %orig value = 0.02

VMLCONF.scal.w_newest = [0.25 0.004 0.0004];
VMLCONF.scal.thres = 0.05;

VMLCONF.plantproj.heights = [250 350 500 700 1000 2000];
VMLCONF.plantproj.extensions = [0 0.05]; %[m] on the 1m plane

VMLCONF.pred.rsun_px = [15 25 45];
VMLCONF.pred.n_cc2pow = 3;
VMLCONF.pred.n_mvec = 3; %1:5
VMLCONF.pred.upd_score = 0.5;
VMLCONF.pred.quantil_score = .1;
VMLCONF.pred.min_nexperts = 0.5;
VMLCONF.pred.tpred = [12:12:300 600 900 1200 1800];
VMLCONF.pred.tpred_persist = 45;

VMLCONF.prj_path = 'E:\ms_proj\';
run([VMLCONF.prj_path 'local_conf.m']);

    