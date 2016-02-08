function c = vmlconfig_cavriglia
c = [];
c.pad_minutes = 10;
c.show_dt = 30;
c.img_ext = '_Debevec.jpeg';
c.basefolder = 'E:\ABB\cav\img\'; %'U:\HDRImages_Cavriglia\img\'; 
%c.basefolder = 'D:\enc_nobak\data\Cavriglia\img\';
c.tmpfolder = 'E:\data\tmp\'; %'D:\enc_nobak\data\Cavriglia\tmp\';
c.datafolder = 'E:\ms_proj\data\Cavriglia\'; %'X:\proj\COSMOS\CloudTracking\Cavriglia\data\';
c.img_offset = [270 250];
c.sky_mask = 'cavriglia_skymask.png';
c.calibration = {'Cavriglia_InternalCalib_20150417.mat',...
  'Cavriglia_Ecalib_20150717.mat',...
  'Cavriglia_model3D.mat'};
c.TimeZone = 'Europe/Zurich';
c.plant_coords = [50,-135,0;195,-135,0;340,-135,0;342,-58,0;345,20,0;188,25,0;32,30,0;40,-53,0]';

c.rsun_px = 150;
c.rsun_px_close = 50;
c.minangle = 5*pi/180;
c.daylight_thres = 0.1;

c.mfi.sz = 300;

c.sundetect.blur_px = 77;
c.sundetect.blur2_px = 17;
c.sundetect.r_roi_px = 50;
c.sundetect.r1_px = 15;
c.sundetect.r2_px = 20;
c.sundetect.star_r1_px = 10;
c.sundetect.star_r2_px = 100;
c.sundetect.satsun_thres = 50;
c.sundetect.star_thres_cov = 0.05;
c.sundetect.star_thres_clear = 0.08;
c.sundetect.maxerr_px = 12;

c.mot.w_cc_newest = 1/3; %weight of the most recent cloud coverage
c.oflow.nframes = 5; %number of recent frames to use
c.oflow.vmax = .10; %max speed [m/sec on the 1m plane]
c.oflow.vstep = 0.005; %initial step [m/sec on the 1m plane]
c.oflow.xtol_rel = 0.01; %relative tolerance for the optimization
c.oflow.heigh_weight_width = 2;   %width of optical flow high weight area [m on plane of height 1]
c.oflow.heigh_weight_cover_center = 1;   %optical flow high weight area needs to cover at least center
c.oflow.low_weight = 0.2; %low weight value (heigh weight value is 1)
c.oflow.t_noise_lim = 1; %limit on noise of time vector
c.sdf.wlo = .01;
c.sdf.vsat = .5;
c.sdf.ub = 10;
c.sdf.sz = 60;
c.sdf.nframes = 10;
c.sdf.alpha = 0.6;

c.seg.ngrid = 15;
c.seg.tile_max_nan = 0.5;
c.seg.alpha = 0.1;
c.seg.minmax_quantil = 0.05;
c.seg.outside_quantil = 0.1;
c.seg.min_peak_mass = 0.15;
c.seg.r2b_resolution = 0.2;
c.seg.r2b_hist_resolution = 0.05;
c.seg.r2b_pad = 0.5;
c.seg.r2b_lpenalty_lim = 0; %no left "outside penalty" right of this value
c.seg.r2b_rpenalty_lim = 0; %no right "outside penalty" left of this value
c.seg.r2b_prior = [-0.5 1];
c.seg.penalty_prior = 0*0.1; %presently disabled
c.seg.v_quadprior = 0;
c.seg.w_quadprior = 0.01;
c.seg.w_local_min = 0.1;
c.seg.w_boundary = 8;
c.seg.r2b_maxrange_singlepeak = 1;
c.seg.wsmooth = 0.2;
c.seg.r2b_dmax = 0.2;
c.seg.rsun_mfi4outer = 80;
c.seg.r2b_multistart_step = .5;
c.seg.max_sunglare_remove = 0.02;
c.seg.locvarpen.thres = [2 3];
c.seg.locvarpen.ksize = 19;
c.seg.locvarpen.maxrange4sky = 0.6;
c.seg.locvarpen.wsky = 0.5;
c.seg.locvarpen.wcloud = 1;

c.scale.default_p_clear = 1e6;
c.scale.default_irr_cloudy = 0.25;
c.scale.warn.p_clear = [.7e6 1.3e6];
c.scale.warn.irr_clear = [.9 1.3];
c.scale.warn.irr_cloudy = [.1 .6];
c.scale.upd.b_irr = 0.25; 
c.scale.upd.b_p = 0.5; 
c.scale.upd.thres_irrel = 0.25;
c.scale.upd.thres_clear = 0.05;

c.plantproj.heights = [250 350 500 700 1000 2000];
c.plantproj.extensions = [0 0.05]; %[m] on the 1m plane

c.pred.tpred = [180 240]; %30:30:900;
c.pred.tpred_persist = 0;
c.pred.tshift = [0 30 60 120];
c.pred.traj = [300 600 900];
c.pred.rsun_px = [15 25 45];
c.pred.mvec_scale = [.95 1 1.05];
c.pred.n_cc2pow = 3;
c.pred.nn_mvec = 3; %1:5
c.pred.eta_p = [0. 0.];
c.pred.eta_irr = [0. 0.];
c.pred.eta_p_agg = [0.1 0.1];
c.pred.eta_irr_agg = [0.1 0.1];
