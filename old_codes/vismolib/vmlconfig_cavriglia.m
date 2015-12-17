VMLCONF = [];
VMLCONF.pad_minutes = 10;
VMLCONF.show_dt = 30;
VMLCONF.img_ext = '.jpeg';
VMLCONF.basefolder = 'U:\HDRImages_Cavriglia\img\'; %'C:\data\cav\img\';
VMLCONF.tmpfolder = 'C:\data\cav\log_data\';
VMLCONF.datafolder = 'C:\Users\chamsei\Documents\GitHub\ms_proj\data\Cavriglia\';
VMLCONF.img_offset = [270 250];
%VMLCONF.sky_mask = 'image_mask_2015_01_21.mat';
VMLCONF.sky_mask = 'cavriglia_skymask.png';
VMLCONF.calibration = {'Cavriglia_InternalCalib_20150417.mat',...
  'Cavriglia_Ecalib_20150717.mat',...
  'Cavriglia_model3D.mat'};
VMLCONF.circumsolarSize = 50;

%% Amrollah conf
VMLCONF.sun_detection_thresold = 15;
VMLCONF.timezone = 'Europe/Rome';
VMLCONF.adapt_window_size = 60; % past seconds which we adapt our irradiation based on
VMLCONF.clear_sky_window = 60*3; % past and future seconds which we look into for clear sky test
VMLCONF.irr_threshold = 5;    % threshold for irradiation  variation to determine clear sky
VMLCONF.irr_comparison_count = 1; % number of elements that we compare in irr log to determine clear sky
VMLCONF.diff_to_reference_irr_threshold = .20;  % thresold to detect cloudy sky based on clear sky refernce
VMLCONF.adaptive_clearsky_reference = 1;
VMLCONF.irr_scale = .98;

VMLCONF.plant_projection_height = 2000; % the height for projection of shape of plant onto image
VMLCONF.plant_coords = [50,-135,0;195,-135,0;340,-135,0;342,-58,0;345,20,0;188,25,0;32,30,0;40,-53,0]';
%VMLCONF.plant_coords = [340,-135,0;345,20,0]';
%%

VMLCONF.cbh = 500;
VMLCONF.clear_sky_power_model = 'Pclearsky.mat';
VMLCONF.rsun_px = 40; % Amri decreased this value from 150
VMLCONF.gaussblur4sundect_px = 3; % changed by Amrollah 77; 
VMLCONF.minangle = 5*pi/180;
VMLCONF.r2b_outliers = 0.01;
VMLCONF.nfilter_mvec = 10;
VMLCONF.nfilter_short_mvec = 5;
VMLCONF.evfilter_short_mvec = 4;
VMLCONF.filter_short_maxstd = 0.2;
VMLCONF.oflow.mask_width = 2;   %width for optical flow mask [m on plane of height 1]
VMLCONF.oflow.mask_cover_center = 1;   %optical flow mask needs to cover at least center
VMLCONF.oflow.masked = 1;       %for optical flow, use image part where clouds come from
VMLCONF.oflow.ignore_diff = 12; %threshold below which red channel differences are ignored
VMLCONF.oflow.w_upwind = 3;     %additional weight for upwind area
VMLCONF.oflow.median_ksize = 3; %KSize for median filter

VMLCONF.thres.nfilter = 5;
VMLCONF.thres.ngrid = 15;
VMLCONF.thres.tile_max_nan = 0.5;
VMLCONF.thres.alpha = 0.1;
VMLCONF.thres.outside_mass = 0.03;  %Amrollah changess it from 0.1
VMLCONF.thres.n_relax_outside = 3;
VMLCONF.thres.r2b_midrange = [0 0.5];
VMLCONF.thres.r2b_resolution = 0.3;
VMLCONF.thres.r2b_hist_resolution = 0.05;
VMLCONF.thres.wsmooth = .1;
VMLCONF.thres.w_local_min = .2;
VMLCONF.thres.w_boundary = 8;
VMLCONF.thres.r2b_dmax = 0.2;
VMLCONF.thres.rsun = 80;
VMLCONF.thres.r2b_multistart_step = .3;



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
