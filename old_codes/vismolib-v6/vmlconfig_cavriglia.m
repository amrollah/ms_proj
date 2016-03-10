function c = vmlconfig_cavriglia
c = [];
c.pad_minutes = 10;
c.show_dt = 30;
c.img_ext = '_Debevec.jpeg';
%c.basefolder = 'U:\HDRImages_Cavriglia\img\';
c.basefolder = 'D:\enc_nobak\data\Cavriglia\img\';
c.tmpfolder = 'D:\enc_nobak\data\Cavriglia\tmp\';
c.datafolder = 'X:\proj\COSMOS\CloudTracking\Cavriglia\data\';
c.img_offset = [270 250];
c.sky_mask = 'cavriglia_skymask.png';
c.calibration = {'Cavriglia_InternalCalib_20150417.mat',...
 'Cavriglia_Ecalib_20150717.mat',...
  'Cavriglia_model3D.mat'};
c.TimeZone = 'Europe/Zurich';
c.plant_coords = [50,-135,0;195,-135,0;340,-135,0;342,-58,0;345,20,0;188,25,0;32,30,0;40,-53,0]';

c.rsun_px = 150;
c.r0sun_px = 40;
c.minangle = 5*pi/180;
c.daylight_thres = 0.1;

c.mfi.sz = 300;

c.sundetect.blur_px = 77;
c.sundetect.blur2_px = 17;
c.sundetect.r_roi_px = 50;
c.sundetect.r0_px = 4;
c.sundetect.r1_px = 15;
c.sundetect.r2_px = 30;
c.sundetect.star_r1_px = 10;
c.sundetect.star_r2min_px = 50;
c.sundetect.star_r2_px = 100;
c.sundetect.blacklevelthres = 50;
c.sundetect.blacklevelthres2 = 150;
c.sundetect.minarel0 = 0.5;
c.sundetect.maxcrit4cov =   [255 1    0.02 0 1];
c.sundetect.mincrit4clear = [150 0.75 0.15 1 0.08];
c.sundetect.star_thres = [0.02 0.05 0.08];
c.sundetect.maxerr_px = 20;

c.corr_sun_bias.clear_frac = 0.1;
c.corr_sun_bias.nfilter = 5;

c.mot.w_cc_newest = 1/3; %weight of the most recent cloud coverage
c.oflow.nframes = 5; %number of recent frames to use
c.oflow.vmax = .10; %max speed [m/sec on the 1m plane]
c.oflow.vstep = 0.005; %initial step [m/sec on the 1m plane]
c.oflow.xtol_rel = 0.01; %relative tolerance for the optimization
c.oflow.heigh_weight_width = 2;   %width of optical flow high weight area [m on plane of height 1]
c.oflow.heigh_weight_cover_center = 1;   %optical flow high weight area needs to cover at least center
c.oflow.low_weight = 0.2; %low weight value (heigh weight value is 1)
c.oflow.t_noise_lim = 1; %limit on noise of time vector
c.sdf.wlo = .001;
c.sdf.vsat = .5;
c.sdf.ub = 10;
c.sdf.sz = 60;
c.sdf.nframes = 11;
c.sdf.nframes_min = 5;
c.sdf.alpha = 0.6;

c.seg.ngrid = 10;
c.seg.rsun_fillclear = 30;
c.seg.tile_max_nan = 0.9;
c.seg.scale_hist = 3;
c.seg.scale_single_peak_crit = 1;
c.seg.single_peak_mse_range = [0.05 0.1];
c.seg.alpha = 0.1;
c.seg.prior_thres = [-1 -0.5 1 1.5];
c.seg.prior_impact = [1 1];
c.seg.minmax_quantil = 0.02;
c.seg.rneigh = 5;
c.seg.r2b_resolution = 0.2;
c.seg.r2b_hist_resolution = 0.05;
c.seg.halfgap = 3;
c.seg.r2b_pad = 0.5;
c.seg.w_boundary = 20;
c.seg.wsmooth = .2;
c.seg.r2b_dmax = 0.5;
c.seg.rsun_maxouter = 400;
c.seg.rsun_close = 180;
c.seg.r2b_multistart_step = .6;
c.seg.max_sunglare_remove = 0.02;
c.seg.locvarthres = [log(4) log(20)];
c.seg.locvarimpact = [0.0 0.2];
c.seg.segcrit_close2sun = [0.01 0.2];
c.seg.rsun_fillclear_precc = 100;

c.scale.default.p_clear = 1e6;
c.scale.default.irr_clear = 1;
c.scale.default.irr_cloudy = 0.25;
c.scale.buflen_min = 3;
c.scale.buflen_max = 9;
c.scale.irrelmin4upd = 0.25;
c.scale.plantccmax4upd = 0.05;

c.plantproj.heights = [250 350 500 700 1000 2000];
c.plantproj.extensions = [0 0.05]; %[m] on the 1m plane

c.pred.tpred = [240]; %30:30:900;
c.pred.tshift = [0 30 60 120];
c.pred.traj = [300 600 900];
c.pred.rsun_px = [20 50];
c.pred.mvec_scale = [.9 1 1.1];
c.pred.t_shapechange = [0 60 120];
c.pred.f_shapechange = [0 0.5 1];
c.pred.cc2class_cutoff = [0.2 0.5 0.8];
c.pred.nn_mvec = 3; %1:5
c.pred.p.sc2w = [240 1];
c.pred.p.sc2w_agg = [240 1];
c.pred.irr.sc2w = [240 1];
c.pred.irr.sc2w_agg = [240 1];
c=local_conf(c);