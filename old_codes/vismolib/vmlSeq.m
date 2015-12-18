classdef vmlSeq < handle
  properties
    files;          %file names
    folder;         %folder
    ti;             %date-time for images
    P;              %PV plant power
    Irr;            %irradiance measurements
    Temp;           %temperature measurements
    conf;           %configuration struct
    oi;             %original image info
    mfi;            %median filtered image info
    sz;             %size of the original image
    calib;          %camera caliblocatio
    strtitle;       %title, typically day as a string
    dt_day;         %day in date time
    dplane;         %min pixel distance on plane
    mvec;           %motion vectors on plane
    printlevel = 1; %0=silent, 1=summarizing, 2=info, 3=verbose
    open_figs;      %list of open figures
    id;             %identification number
    varname;        %variable name in base workspace
    Z_tpred         %prediction times that feature vectors are available for
    Z;              %cell array of feature vectors
    Z_npoints;      %cell array of avg / min number of points used to compute each Z entry
    W_PBP;          %pairs prediction times / blurs where trained weights are available 
    W;              %cell array of trained weights
    W_RMSE;         %training RMSE
    Pclearsky;      %clear sky power model
    thres;          %data for cloud thresholding
    xcur;           %current frame for feature computation
    ext_calib;      %external calibration matrix (R)
    sun_pos;        % adjusted sun positions for all images
    ClearSkyRef;    % PVL reference clear-sky global,direct,diffuse irradiation
    ClearSkyOrigRef;  % PVL reference clear-sky model without adaptation
    is_clear_states; % clearness state for all times
    clear_times;    % index of clear times in ti matrix
    cloudy_times;   % index of cloudy times in ti matrix
    difff;
    adjuster_ei;
  end
  
  methods
  
%% basic functions: constructor, imread, getP, ...

    function obj = vmlSeq(folder,hours,conf)
      if nargin<2, hours = []; end
      if nargin<3, conf = evalin('base','VMLCONF'); end
      
      %each vmlSeq object gets an id, in order to identify it back
      %from Matlab figures
      obj.id = floor(rand*1e12);
      
      %store configuration and list image files
      obj.conf = conf;
      [obj.files, obj.folder, obj.ti] = vmlListImgFiles(folder, hours, conf);
      
      %read in sky mask
      if isempty(conf.sky_mask)
        x1 = obj.imread(1);
        obj.oi.sm = true(size(x1,1),size(x1,2));
      elseif strfind(conf.sky_mask,'.png')
        x1 = imread([conf.datafolder conf.sky_mask]);
        obj.oi.sm = (x1(:,:,1)>0);
      else
        x1=load([conf.datafolder conf.sky_mask]);
        fn = fieldnames(x1);
        obj.oi.sm = x1.(fn{1});
        obj.oi.sm = obj.oi.sm(:,:,1);
      end
      
      %read in internal / external calibration and location info
      obj.calib = [];
      x1=load([conf.datafolder conf.calibration{1}]);
      fn = fieldnames(x1);
      obj.calib.intern = x1.(fn{1});
      if isempty(conf.calibration{2}), obj.calib.Rext = eye(3);
      else
        x1=load([conf.datafolder conf.calibration{2}]);
        fn = fieldnames(x1);
        obj.calib.Rext = x1.(fn{1});
      end
      x1=load([conf.datafolder conf.calibration{3}]);
      fn = fieldnames(x1);
      obj.calib.model3D = x1.(fn{1});
            
      %compute pixels mapped to unit sphere and indentify / cut away
      %pixels too close to the horizon
      [height,width] = size(obj.oi.sm);
      [yy,xx] = ndgrid(1:height,1:width);
      obj.oi.spx = obj.im2camworld([yy(:)';xx(:)'],[0 0]);
      obj.oi.sm(obj.oi.spx(3,:)<sin(conf.minangle)) = false;
      obj.oi.spx = reshape(obj.oi.spx',[height width 3]);

      %now the final area of interest (bounding box) can be determined
      obj.oi.sky_area = [max(0,find(any(obj.oi.sm,2),1)-1) ...
        min(height,find(any(obj.oi.sm,2),1,'last')+1) ...
        max(0,find(any(obj.oi.sm,1),1)-1) ...
        min(width,find(any(obj.oi.sm,1),1,'last')+1)];

      %reduce sky mask to final bounding box
      obj.oi.sm = obj.oi.sm(obj.oi.sky_area(1):obj.oi.sky_area(2),obj.oi.sky_area(3):obj.oi.sky_area(4));
      obj.sz = size(obj.oi.sm);
      [obj.oi.YY,obj.oi.XX] = ndgrid(1:obj.sz(1),1:obj.sz(2));
      
      %compute border of the sky mask
      obj.oi.sm_border = false(obj.sz);
      obj.oi.sm_border(1:end-1,:) = obj.oi.sm_border(1:end-1,:) | obj.oi.sm(2:end,:);
      obj.oi.sm_border(2:end,:) = obj.oi.sm_border(2:end,:) | obj.oi.sm(1:end-1,:);
      obj.oi.sm_border(:,1:end-1) = obj.oi.sm_border(:,1:end-1) | obj.oi.sm(:,2:end);
      obj.oi.sm_border(:,2:end) = obj.oi.sm_border(:,2:end) | obj.oi.sm(:,1:end-1);
      obj.oi.sm_border(obj.oi.sm) = false;
      
      %reduce pixels on unit sphere to final bounding box and map
      %them to unit plane
      obj.oi.spx = obj.oi.spx(obj.oi.sky_area(1):obj.oi.sky_area(2),obj.oi.sky_area(3):obj.oi.sky_area(4),:);
      obj.oi.spx = reshape(obj.oi.spx,numel(obj.oi.sm),3);
      obj.oi.ppx = obj.oi.spx(:,1:2)./repmat(obj.oi.spx(:,3),[1 2]);
      obj.oi.ppx(~(obj.oi.sm | obj.oi.sm_border),:) = inf;
      
      %compute minimum distances of neighboring pixels on unit plane
      plane_px1 = reshape(obj.oi.ppx,[obj.sz 2]);
      d1 = plane_px1(2:end,:,1)-plane_px1(1:end-1,:,1);
      d2 = plane_px1(:,2:end,2)-plane_px1(:,1:end-1,2);
      obj.dplane = min(min(abs(d1(~isnan(d1)))),min(abs(d2(~isnan(d2)))));
      assert(obj.dplane >= 1/max(height,width));
      
      
      
      %load irradiation and power data
      data_Power = [];
      try
        run([conf.datafolder 'loaddata.m']);
      catch
          disp('Error in loading log files.');
      end
      crop_fun = @(x)x(max(1,find(x(:,1)>=min(obj.ti),1)-1):min(size(x,1),...
        find(x(:,1)<=max(obj.ti)+conf.pad_minutes*60/86400,1,'last')+1),:);
      if ~isempty(data_Power)
          obj.P = crop_fun(data_Power); 
          P_tmp = sort(obj.P(:,2));
          if mean(P_tmp(1:500))< -1000
           disp('Invalid power log files.');
          end
      end 
      if ~isempty(data_Irr), obj.Irr = crop_fun(data_Irr); end
      if ~isempty(data_Temp), obj.Temp = crop_fun(data_Temp); end
      
      %assign title and day in date time format
      obj.strtitle = folder;
      obj.dt_day = datenum(obj.strtitle,'yyyy_mm_dd');
         
      %initialize data for motion vectors and power prediction
      obj.mvec = nan(2,length(obj.ti));
      obj.Z_tpred = [];
      obj.Z = {};
      obj.Z_npoints = {};
      obj.W_PBP = [];
      obj.W = {};
      obj.W_RMSE = [];
      
      %load clear sky power model
      if ~isempty(conf.clear_sky_power_model)
        x1=load([conf.datafolder conf.clear_sky_power_model]);
        fn = fieldnames(x1);
        x1 = x1.(fn{1});
        %day = obj.dt_day-datenum(obj.strtitle(1:4),'yyyy');
        obj.Pclearsky = x1;
      end
      
      if ~isempty(obj.Irr)
          LinkeTurbidity = prep_LinkeTurbidity(obj);
          UTC_offset = obj.calib.model3D.UTC;
          if ~isdst(datetime(datevec(obj.ti(end)),'TimeZone',obj.conf.timezone))
            UTC_offset =  UTC_offset - 1;
          end
          times = pvl_maketimestruct(obj.Irr(:,1), UTC_offset);
          [ClearSkyGHI,ClearSkyDNI,ClearSkyDHI] = pvl_clearsky_ineichen(times,obj.calib.model3D.Location,'LinkeTurbidityInput',LinkeTurbidity);

          obj.ClearSkyOrigRef = [obj.Irr(:,1),ClearSkyGHI*obj.conf.irr_scale,ClearSkyDNI*obj.conf.irr_scale,ClearSkyDHI*obj.conf.irr_scale]; 
          obj.ClearSkyRef = obj.ClearSkyOrigRef;
          % amrollah: clear-sky detection and cloudy-sky detection based on irradiation
          % log file
          obj.is_clear_states = arrayfun(@(t) is_clear(obj, t), obj.Irr(:,1));
          obj.clear_times = find(obj.is_clear_states==1);
          obj.cloudy_times = find(obj.is_clear_states<1);
          
          % adaptive irradiation reference
          %% TEMP Comment
          obj.ClearSkyRef = [obj.Irr(:,1),adjust_reference(obj,ClearSkyGHI),ClearSkyDNI.*obj.conf.irr_scale,ClearSkyDHI.*obj.conf.irr_scale];
      end
      
      obj.ext_calib = load([conf.datafolder conf.calibration{2}]);
      % calculate sun positions for all images and store them in object
      %% TEMP Comment
      obj.calc_sun_positions();
      
      %do a median downscale of the image and compute the sky mask...
      [obj.mfi.sm,obj.mfi.xx,obj.mfi.yy] = ...
        vmlMedianDownscale(obj.oi.sm,obj.conf.thres.nfilter);
      obj.mfi.sm = logical(obj.mfi.sm);
      
      %... and the unit plane pixels for the median downscaled images 
      [obj.mfi.YY,obj.mfi.XX] = ndgrid(obj.mfi.yy,obj.mfi.xx);
      spx = obj.im2camworld([obj.mfi.YY(:)';obj.mfi.XX(:)'])';
      obj.mfi.ppx = spx(:,1:2)./repmat(spx(:,3),[1 2]);
      obj.mfi.ppx(~obj.mfi.sm,:) = nan;
      
     %tiles for the local thresholding
      [sy,sx] = size(obj.mfi.sm);
      rr = ((1:obj.conf.thres.ngrid)-.5)/obj.conf.thres.ngrid;
      obj.thres.xxg = [1 round(1+(sx-1)*rr) sx]; 
      obj.thres.yyg = [1 round(1+(sy-1)*rr) sy]; 
      obj.thres.have_v = zeros(obj.conf.thres.ngrid);
      for ix=1:obj.conf.thres.ngrid
        for iy=1:obj.conf.thres.ngrid
          if any(any(obj.mfi.sm(obj.thres.yyg(iy):obj.thres.yyg(iy+2),obj.thres.xxg(ix):obj.thres.xxg(ix+2))))
            obj.thres.have_v(iy,ix) = 1;
          end
        end
      end
      rr = (0:obj.conf.thres.ngrid)/obj.conf.thres.ngrid;
      obj.thres.xxgb = round(sx*rr); 
      obj.thres.yygb = round(sy*rr); 
      yy = obj.thres.yyg(2:end-1)-(sy+1)/2;
      xx = obj.thres.xxg(2:end-1)-(sx+1)/2;
      [YY,XX]=ndgrid(yy,xx);
      obj.thres.isouter = XX.^2+YY.^2>=min(max(abs(yy)),max(abs(xx)))^2;
    end
    
    function x = imread(obj,j)
      %read a single frame
      x = imread([obj.folder obj.files{j} obj.conf.img_ext]);
      if ~isempty(obj.oi.sky_area)
        x = x(obj.oi.sky_area(1):obj.oi.sky_area(2),obj.oi.sky_area(3):obj.oi.sky_area(4),:);
      end
    end
        
    function P1 = getP(obj,t)
      %get the power value for time(s) t
      P1 = interp1(obj.P(:,1),obj.P(:,2),t);
    end
    
    function P1 = getPclearsky(obj,t)
      %compute the clear sky power model for time(s) t
      tt = t(:)'-floor(obj.ti(1));
      P1 = 1;%sim_nn_nf(obj.Pclearsky.w,tt);
    end
    
    function [jZ, tpred] = getZidx(obj,tpred,accept_empty)
      if nargin<3, accept_empty = 0; end
      if isempty(obj.Z_tpred)
        if accept_empty, jZ = []; return;
        else error('no features computed yet'); end
      end
      if isempty(tpred), tpred = obj.Z_tpred(1); end
      jZ = find(tpred==obj.Z_tpred);
      if isempty(jZ) && ~accept_empty, error('no features computed for this tpred'); end
    end
    
    function [jW, tpred, tblur, usepow] = getWidx(obj,PBP,accept_empty)
      if nargin<3, accept_empty = 0; end
      if isempty(obj.Z_tpred), error('no features computed yet'); end
      if isempty(obj.W_PBP)
        if accept_empty, jW = []; return;
        else error('no weights trained yet'); end
      end
      if isempty(PBP), tpred = obj.Z_tpred(1); tblur = [];
      else tpred = PBP(1); tblur = PBP(2); usepow = PBP(3); end
      if isempty(tblur), jW = find(obj.W_PBP(1,:)==tpred,1);
      else jW = find(obj.W_PBP(1,:)==tpred & obj.W_PBP(2,:)==tblur & obj.W_PBP(3,:)==usepow); end
      if isempty(jW) && ~accept_empty, error('no weights for this PBP'); end
      if ~isempty(jW),  tblur = obj.W_PBP(2,jW); usepow = obj.W_PBP(3,jW); end
    end
    
    function data = export(obj)
      data = [];
      data.ti = obj.ti;
      data.mvec = obj.mvec;
      data.Z_tpred = obj.Z_tpred;
      data.Z = obj.Z;
      data.Z_npoints = obj.Z_npoints;
    end
    
    function import(obj,data)
      [~,j,j1] = intersect(obj.ti,data.ti);
      obj.mvec(:,j) = data.mvec(:,j1);
      obj.Z_tpred = data.Z_tpred;
      N = length(obj.Z_tpred);
      obj.Z = cell(1,N);
      obj.Z_npoints = cell(1,N);
      for i = 1:N
        obj.Z{i} = nan(obj.conf.pred.n_feature_maps*obj.conf.pred.ngrid^2,length(obj.ti));
        obj.Z_npoints{i} = nan(2,length(obj.ti));
        if size(obj.Z{i},1)==size(data.Z{i},1)
          obj.Z{i}(:,j) = data.Z{i}(:,j1);
          obj.Z_npoints{i}(:,j) = data.Z_npoints{i}(:,j1);
        end
      end
    end
    
    function j = getidx(obj,t1,t2)
      %get frame number(s) for given time(s)
      if nargin<3, t2=[]; end
      if ischar(t1), t1 = time_str2num(t1); end
      if ischar(t2), t2 = time_str2num(t2); end
      if isempty(t2), [~,j] = min(abs(obj.ti-floor(obj.ti(1))-t1/24));
      else j = find(obj.ti-floor(obj.ti(1))>=t1/24 & obj.ti-floor(obj.ti(1))<t2/24); end
    end
    
    function l = cloudiness_label(obj, tt)
      l = zeros(1,length(tt));
      ff = obj.getP(tt)./obj.getPclearsky(tt);
      l(ff>=obj.conf.pred.clearsky_thres) = 1;
      l(ff<=obj.conf.pred.covered_thres) = -1;
    end
    
    function loadframe(obj,j)
      %load frame number j:
      % - read the image
      % - compute the median downscaled image
      % - identify opposite class neighbors
      % - optimize the cloud segmentation threshold map
      if j<1 || j>length(obj.ti), obj.xcur=[]; return; end
      if ~isempty(obj.xcur) && any(obj.xcur.j == j), return; end
      obj.xcur.j = [];
      obj.xcur.x0 = obj.imread(j);
      obj.xcur.x = vmlMedianDownscale(obj.xcur.x0,obj.conf.thres.nfilter);
      obj.xcur.r2b0 = vmlRed2Blue(obj.xcur.x0);
      x1 = obj.xcur.r2b0(obj.oi.sm);
      obj.xcur.r2b = vmlMedianDownscale(obj.xcur.r2b0,obj.conf.thres.nfilter);
      obj.xcur.lum = vmlMedianDownscale(mean(obj.xcur.x0,3)/255,obj.conf.thres.nfilter);
      
      yx = sunpos_im(obj,j);
      b = ((obj.mfi.XX-yx(2)).^2+(obj.mfi.YY-yx(1)).^2<=obj.conf.rsun_px^2);
      obj.xcur.avglum = mean(obj.xcur.lum(logical(obj.mfi.sm)));
      obj.xcur.avglumsun = mean(obj.xcur.lum(b));
      obj.xcur.sm = obj.mfi.sm & ~b;
      
      r2b_midrange = obj.conf.thres.r2b_midrange;
      r2b = obj.xcur.r2b; r2b(~obj.xcur.sm) = NaN;
      r2bmin = min(r2b_midrange(1),min(obj.xcur.r2b(obj.mfi.sm)));
      r2bmax = max(r2b_midrange(2),max(obj.xcur.r2b(obj.mfi.sm)));
      res = obj.conf.thres.r2b_resolution;
      resh = obj.conf.thres.r2b_hist_resolution;
      N = round((r2bmax-r2bmin)/res)+2;
      obj.xcur.thres.zz = r2bmin+(r2bmax-r2bmin-N*res)/2+res*(0:N);
      obj.xcur.thres.zzh = r2bmin-resh*1.5:resh:r2bmax+resh*2;
      NN = length(obj.xcur.thres.zzh);
      jmid_lo = find(obj.xcur.thres.zzh>=r2b_midrange(1),1);
      jmid_hi = find(obj.xcur.thres.zzh<=r2b_midrange(2),1,'last');
      obj.xcur.thres.X = tanh(bsxfun(@minus,obj.xcur.thres.zzh',obj.xcur.thres.zz)/res);
      obj.xcur.thres.ww0 = nan(N+1,obj.conf.thres.ngrid,obj.conf.thres.ngrid);
      obj.xcur.thres.ww = nan(N+1,obj.conf.thres.ngrid,obj.conf.thres.ngrid);
      obj.xcur.thres.bb = nan(obj.conf.thres.ngrid);
      obj.xcur.thres.j0 = nan(obj.conf.thres.ngrid);
      obj.xcur.thres.isouter = obj.thres.isouter;
      obj.xcur.thres.histn = zeros(NN,obj.conf.thres.ngrid,obj.conf.thres.ngrid);
      obj.xcur.thres.v0_mask = zeros(size(obj.xcur.sm));
      sunp = obj.sunpos_im(j);
      [~,sunpy]=min(abs(obj.mfi.yy-sunp(1))); 
      [~,sunpx]=min(abs(obj.mfi.xx-sunp(2)));
      sunp = [sunpy;sunpx];
      for ix=1:obj.conf.thres.ngrid
        for iy=1:obj.conf.thres.ngrid
          if sum(([obj.thres.yyg(iy+1);obj.thres.xxg(ix+1)]-sunp).^2) <= ...
              obj.conf.thres.rsun^2
            obj.xcur.thres.isouter(iy,ix) = 0;
          end
          z = r2b(obj.thres.yyg(iy):obj.thres.yyg(iy+2),obj.thres.xxg(ix):obj.thres.xxg(ix+2));
          z = z(:);
          jjy = obj.thres.yygb(iy)+1:obj.thres.yygb(iy+1); 
          jjx = obj.thres.xxgb(ix)+1:obj.thres.xxgb(ix+1);
          if mean(isnan(z))>obj.conf.thres.tile_max_nan
            obj.xcur.thres.v0_mask(jjy,jjx) = -2;
          else
            z = z(~isnan(z));
            nn=hist(z,obj.xcur.thres.zzh)';
            w = wlsq(obj.xcur.thres.X,nn,[],...%.1+sqrt(mean(nn)./max(nn,1)),...
              obj.conf.thres.alpha*ones(1,N+1),[],[],...
              -obj.xcur.thres.X,nn*0,ones(1,N+1),0);
            yy = obj.xcur.thres.X*w;
            scale = sum(nn)/prod(2*size(obj.mfi.sm)/obj.conf.thres.ngrid)/sum(yy);
            yy = yy*scale; w = w*scale; nn = nn*scale;
            j0 = otsu_local_min(yy);
            if min(sum(yy(1:j0-1)),sum(yy(j0+1:end)))<sum(yy)*obj.conf.thres.outside_mass
              j0 = NaN;
            end
            obj.xcur.thres.j0(iy,ix) = j0;
            obj.xcur.thres.ww0(:,iy,ix) = w;
            obj.xcur.thres.histn(:,iy,ix) = nn;
            if isnan(j0)
              obj.xcur.thres.v0_mask(jjy,jjx) = -1;
            else
              obj.xcur.thres.v0_mask(jjy,jjx) = obj.xcur.r2b(jjy,jjx)>=obj.xcur.thres.zzh(j0);
            end
            %repeat the fitting, but with high penalties at both ends
            b = max(nn);
            nnavg = mean(nn);
            s1 = sum(nn)*obj.conf.thres.outside_mass;
            j1 = max(1,min(jmid_lo,find(cumsum(nn)>=s1,1))-obj.conf.thres.n_relax_outside); 
            j2 = max(1,min(NN+1-jmid_hi,find(cumsum(nn(end:-1:1))>=s1,1))-obj.conf.thres.n_relax_outside);
            nn(1:j1) = max(nn(1:j1),b+(nnavg-b)*(0:j1-1)'/j1);
            nn(NN:-1:NN+1-j2) = max(nn(NN:-1:NN+1-j2),b+(nnavg-b)*(0:j2-1)'/j2);
            w = -wlsq(obj.xcur.thres.X,b-nn,[],...
              obj.conf.thres.alpha*ones(1,N+1),[],[],...
              [-obj.xcur.thres.X;obj.xcur.thres.X],[nn*0;nn*0+b],ones(1,N+1),0);
            obj.xcur.thres.ww(:,iy,ix) = w;
            obj.xcur.thres.bb(iy,ix) = b;
          end
        end
      end
      
      obj.xcur.thres.v = vmlSmoothThresMap(obj.xcur.thres,obj.thres.have_v,obj.conf.thres);
      [YYg,XXg] = ndgrid(obj.thres.yyg,obj.thres.xxg);
      VVg = XXg;
      VVg(2:end-1,2:end-1) = obj.xcur.thres.v;
      VVg(1,:) = VVg(2,:); VVg(end,:) = VVg(end-1,:);
      VVg(:,1) = VVg(:,2); VVg(:,end) = VVg(:,end-1);
      [YY,XX] = ndgrid(1:length(obj.mfi.yy),1:length(obj.mfi.xx));
      obj.xcur.thres.vmfi = interp2(XXg,YYg,VVg,XX,YY);
      obj.xcur.thres.vmfi(~obj.mfi.sm) = NaN;

      obj.xcur.j = j;      
    end
      
      
% end basic functions

%% plotting functionality
    
    function newfig(obj,figno,tit)
      if isempty(obj.varname)
        obj.varname = '<local>';
        for v = evalin('base','who')' 
          if strcmp(evalin('base',['class(' v{1} ')']),'vmlSeq');
            if (evalin('base',[v{1} '.id'])==obj.id)
              obj.varname = v{1};
              break;
            end
          end
        end
      end
      obj.open_figs = unique([obj.open_figs figno]);
      figure(figno); clf;
      set(gcf,'Name',[tit ' @@@ ' obj.varname]);%,'NumberTitle','off')
    end
    
    function showSun(obj,j,col,marker)
      if nargin<3, col = 'k'; end
      if nargin<4, marker = []; end
      tt = (0:.1:360)*pi/180;
      yx = sunpos_im(obj,j);
      xx = yx(2)+cos(tt)*obj.conf.rsun_px;
      xx(xx<1 | xx>obj.sz(2)) = NaN;
      yy = yx(1)+sin(tt)*obj.conf.rsun_px;
      yy(yy<1 | yy>obj.sz(1)) = NaN;
      if ~isempty(marker), line(yx(2),yx(1),'color',col,'marker',marker); end
      line(xx,yy,'color',col);
    end
    
    function showframe(obj,j,show_sun_or_skymask)
      %show the frame, part which does not belong to the sky is
      %greened out
      if nargin<3, show_sun_or_skymask = 0; end
      if ~isscalar(show_sun_or_skymask), sm = show_sun_or_skymask;
      elseif show_sun_or_skymask>=2, sm = obj.skymask_wo_sun(j);
      else sm = obj.oi.sm; 
      end
      x = vmlColorify(obj.imread(j),~sm,2,64);
      image(x); axis off;
      if isscalar(show_sun_or_skymask) && show_sun_or_skymask==1
        obj.showSun(j,'g','.');
      end
      title([datestr(obj.ti(j),'HH:MM:SS') ' (#' num2str(j) ')']);
    end
    
    function show_ecalib(obj)
      %produces plot to visually check the external calibration
      if ~isfield(obj.calib,'kabsch') || isempty(obj.calib.kabsch)
        error('no external calibration was executed on this vmlSeq object');
      end
      obj.newfig(1329,'External calibration');
      subplot(3,2,5:6);
      plot(obj.calib.kabsch.dt,obj.calib.kabsch.im_model(1,:),'-b.',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_model(2,:),'-r.',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_detect(1,:),'bo',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_detect(2,:),'ro');
      title(obj.strtitle,'interpreter','none');
      datetickzoom;
      show_ecalib_upd(1);
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@show_ecalib_upd)
      function out = show_ecalib_upd( j, event_obj)
        show_processed = 0;
        if nargin>=2
          [~,j] = min(abs(event_obj.Position(1)-obj.calib.kabsch.dt));
          show_processed = (diff(abs(obj.calib.kabsch.im_model(:,j)-event_obj.Position(2)))<0);
        end
        figure(1329);
        subplot(3,2,1:4); 
        if ~show_processed
          obj.showframe(obj.calib.kabsch.frame_nos(j));
        else
          x = obj.im4sundetect(obj.calib.kabsch.frame_nos(j));
          imagesc(x); axis off;
          title([datestr(obj.ti(obj.calib.kabsch.frame_nos(j)),'HH:MM:SS') ' (#' num2str(obj.calib.kabsch.frame_nos(j)) ')']);
        end
        hold on; 
        plot(obj.calib.kabsch.im_model(2,:),obj.calib.kabsch.im_model(1,:),'-r.',...
          obj.calib.kabsch.im_detect(2,:),obj.calib.kabsch.im_detect(1,:),'ro');
        hold off;
        drawnow;
        if nargin<2, out = [];
        else out = datestr(event_obj.Position(1),'HH:MM:SS'); end
      end
    end
    
    function plotpow(obj, toffset)
      %plot the power profiles
      if nargin<2, toffset=0; end
      toffset = toffset/86400;
      l = obj.cloudiness_label(obj.ti+toffset);
      plot(obj.P(:,1),obj.getP(obj.P(:,1)+toffset),'b',...
        obj.ti,obj.getP(obj.ti+toffset),'b.',...
        obj.ti(l>0),obj.getP(obj.ti(l>0)+toffset),'r.',...
        obj.ti(l<0),obj.getP(obj.ti(l<0)+toffset),'k.',...
        obj.P(:,1),obj.getPclearsky(obj.P(:,1)+toffset),'r');
      datetickzoom;
    end
    
    function plot_class_pred(obj,tpred,class_pred)
      l = obj.cloudiness_label(obj.ti);
      Phi = mean(obj.getP(obj.ti(l>0)));
      Plo = mean(obj.getP(obj.ti(l<0)));
      obj.plotpow(tpred); hold on;
      plot([min(obj.ti) max(obj.ti)],(Phi+Plo)/2+[0 0],'m:',...
        obj.ti,(Phi+Plo)/2+saturate(class_pred)*(Phi-Plo)/2,'m'); 
      hold off;
    end
    
    function showtraj(obj,j,tpred,use_filt)
      %plot the motion trajectory for frame j
      if nargin<3, tpred=0; end
      if nargin<4, use_filt=1; end
      obj.showframe(j,1);
      if use_filt, mvec1 = obj.mvec_filtered(j);
      else mvec1 = obj.mvec(:,j); end
      if all(~isnan(mvec1))
        n = 300; n1 = 60;
        p = obj.sunpos_plane(j);
        pp = repmat(p,1,n)-mvec1*((1:n));
        yx = obj.camworld2im([pp;ones(1,n)]);
        yx_tpred = obj.camworld2im([p-mvec1*tpred;1]);
        hold on;
        plot(yx(2,:),yx(1,:),'r',yx(2,n1:n1:n),yx(1,n1:n1:n),'ro',...
          yx_tpred(2),yx_tpred(1),'rx');
        hold off;
      end
      title([datestr(obj.ti(j),'HH:MM:SS') ' #' num2str(j) ' tpred=' num2str(tpred)]);
    end
    
    function showLum(obj,j,~)
      obj.loadframe(j);
      x = obj.xcur.lum;
      x(~obj.mfi.sm) = NaN;
      imagesc(obj.mfi.yy,obj.mfi.xx,x); axis off;
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5 ]);
      %set(imagesc(x),'alphadata',~isnan(x));axis off;
    end
    
    function showR2B0(obj,j,~)
      obj.loadframe(j);
      x = obj.xcur.r2b;
      x(~obj.mfi.sm) = NaN;
      h = imagesc(obj.mfi.yy,obj.mfi.xx,x); axis off;
      set(h,'alphadata',~isnan(x));
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5 ]);
    end
    
    function showR2B(obj,j,~)
      obj.loadframe(j);
      x = vmlColorify(vmlColorify(repmat(uint8(obj.xcur.r2b*192),[1 1 3]),...
        obj.xcur.cloudy,1,64),~obj.xcur.sm,2,64);
      image(obj.mfi.yy,obj.mfi.xx,x); axis off;
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5 ]);
      tt = (0:.1:360)*pi/180;
      tt1 = (0:obj.conf.thres.nsec-1)/obj.conf.thres.nsec*2*pi;
      xc = obj.calib.intern.yc-obj.oi.sky_area(3)+1;
      yc = obj.calib.intern.xc-obj.oi.sky_area(1)+1;
      hold on; 
      plot(xc+sin(tt)*obj.conf.thres.rinner,yc+cos(tt)*obj.conf.thres.rinner,'g',...
        xc+sin(tt)*obj.conf.thres.router,yc+cos(tt)*obj.conf.thres.router,'g',...
        xc+sin(tt1)*obj.conf.thres.router,yc-cos(tt1)*obj.conf.thres.router,'go'); 
      hold off;
    end
    
    function showNeigh(obj,j,~)
      if nargin<2, j = obj.xcur.j; end
      obj.loadframe(j);
      x = obj.xcur.x;
      m = false(size(obj.xcur.r2b));
      m(obj.xcur.neighInfo) = true;
      x = vmlColorify(x,m,2,32);
      m = false(size(obj.xcur.r2b));
      m(obj.xcur.nneighInfo) = true;
      x = vmlColorify(x,m,1,64);
      image(obj.mfi.yy,obj.mfi.xx,x); axis off;
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5 ]);
    end
    
function showThres(obj,j)
      if nargin<2, j = obj.xcur.j;
      else obj.loadframe(j); end
      obj.newfig(1328,['Cloud thresholds, frame #' num2str(j) ' / ' datestr(obj.ti(j))]);
      cb_text = {'r2b','r2b col','thres0','thres','thres val','contour'};
      hcb = zeros(1,length(cb_text));
      htxt = uicontrol('Style','text','Position',[10 5 650 16]);
      for i = 1:length(hcb)
        hcb(i) = uicontrol('Style','checkbox','Value',i==length(hcb),...
          'Position',[i*85 24 80 16],'String',cb_text{i});
      end
      auxplot_shown = [-1 -1];
      set(gcf,'toolbar','figure');
      set(hcb,'Callback',{@showThres_upd}); %,hcb,hax
      showThres_upd;

      function showThres_click(varargin)
        p = get(gca,'currentpoint');
        [~,jx] = min(abs(p(1,1)-obj.mfi.xx));
        [~,jy] = min(abs(p(1,2)-obj.mfi.yy));
        ix = max(1,min(obj.conf.thres.ngrid,1+floor((p(1,1)-1)/(obj.sz(2)-1)*obj.conf.thres.ngrid)));
        iy = max(1,min(obj.conf.thres.ngrid,1+floor((p(1,2)-1)/(obj.sz(1)-1)*obj.conf.thres.ngrid)));
        set(htxt,'String',['xy = [' num2str(obj.mfi.xx(jx)) ', ' num2str(obj.mfi.yy(jy)) ...
           '] RGB = (' strjoin(cellfun(@num2str,squeeze(num2cell(obj.xcur.x(jy,jx,:))),'UniformOutput',0)',', ') ...
           ') R2B = ' num2str(obj.xcur.r2b(jy,jx),'%.2f') ' / ' ...
           'xytile = [' num2str(ix) ', ' num2str(iy) ']']);
         if any([ix iy]~=auxplot_shown)
           nn = obj.xcur.thres.histn(:,iy,ix);
           zzh = obj.xcur.thres.zzh;
           yy0 = obj.xcur.thres.X*obj.xcur.thres.ww0(:,iy,ix);
           yy = obj.xcur.thres.X*obj.xcur.thres.ww(:,iy,ix)+obj.xcur.thres.bb(iy,ix);
           j0 = obj.xcur.thres.j0(iy,ix);
           figure(13281); clf; hold on;
           bar(zzh,nn,1); 
           plot(zzh,max(0,yy0),'r',zzh,max(0,yy),'g--','linewidth',2);
           if ~isnan(j0), plot(zzh(j0),max(0,yy0(j0)),'ro','linewidth',2); end
           plot(obj.xcur.thres.v(iy,ix),0,'gd','linewidth',2);
           hold off; grid on;
         end
      end
      
      function showThres_upd(varargin)
        figure(1328);
        vv = get(hcb,'value');
        if vv{2} 
          x = obj.xcur.r2b;
          x(~obj.mfi.sm) = 0;
          ih = imagesc(obj.mfi.xx,obj.mfi.yy,x); axis off;
          pos = get(gca,'Position');
          colorbar('vert');
          set(gca,'Position',pos);
        else
          if ~vv{1}, x = obj.xcur.x;
          else
            x = obj.xcur.r2b; x = x-min(x(:));
            x = repmat(uint8(x/max(x(:))*192),[1 1 3]); 
          end
          x = vmlColorify(x,~obj.mfi.sm,1:3,-255);
          if vv{4}
            x = vmlColorify(x,obj.xcur.r2b<obj.xcur.thres.vmfi,1,-40); 
            x = vmlColorify(x,obj.xcur.r2b>=obj.xcur.thres.vmfi,1,40); 
          elseif vv{3}
            x = vmlColorify(x,obj.xcur.thres.v0_mask==-2,1:3,-80); 
            x = vmlColorify(x,obj.xcur.thres.v0_mask==-1,1:3,40);
            x = vmlColorify(x,obj.xcur.thres.v0_mask==0,1,-40); 
            x = vmlColorify(x,obj.xcur.thres.v0_mask==1,1,40); 
          end
          ih = image(obj.mfi.xx,obj.mfi.yy,x); axis off;
        end
        title([datestr(obj.ti(j),'HH:MM:SS') ' (#' num2str(j) ')']);
        axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5]);
        set(ih,'ButtonDownFcn',@showThres_click)
        obj.showSun(j);
        for p=(1:obj.conf.thres.ngrid-1)/obj.conf.thres.ngrid
          line(1+[p p]*(obj.sz(2)-1),[1 obj.sz(1)],'Color',[0 .5 0]);
          line([1 obj.sz(2)],1+[p p]*(obj.sz(1)-1),'Color',[0 .5 0]);
        end
        for iy=1:obj.conf.thres.ngrid
          for ix=1:obj.conf.thres.ngrid
            if ~obj.thres.have_v(iy,ix)
              line(1+(ix-.5)/obj.conf.thres.ngrid*(obj.sz(2)-1),1+(iy-.5)/obj.conf.thres.ngrid*(obj.sz(1)-1),'Marker','x','Color',[0 .5 0])
            end
          end
        end
        if vv{5} && ~vv{2}
          hold on;
          [~,ih] = contour(obj.mfi.XX,obj.mfi.YY,obj.xcur.thres.vmfi,20);
          set(gca,'clim',[min(obj.xcur.thres.vmfi(:)) max(obj.xcur.thres.vmfi(:))]);
          hold off;
          set(ih,'ButtonDownFcn',@showThres_click)
          pos = get(gca,'Position');
          colorbar('vert');
          set(gca,'Position',pos);
        elseif vv{6}
          hold on;
          [~,ih] = contour(obj.mfi.XX,obj.mfi.YY,obj.xcur.r2b-obj.xcur.thres.vmfi,[0 0],char('g'+vv{2}*('k'-'g')),'linewidth',1+vv{2});
          hold off;
          set(ih,'ButtonDownFcn',@showThres_click)
        end
        set(gca,'xlim',[0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5]);
        set(gca,'ylim',[0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5]);
      end
    end
        
    function showThresVmap(obj,j,~)
      obj.loadframe(j);
      contour(obj.mfi.XX,-obj.mfi.YY,obj.xcur.thresvmap,60); 
      hold on; 
      contour(obj.mfi.XX,-obj.mfi.YY,obj.xcur.thresvmap,obj.xcur.thresv0+[0 0],'k');
      hold off;
      set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
      title(['#' num2str(j) ' ' num2str(min(obj.xcur.thresvmap(:)),'%.2f') ' - ' ...
        num2str(max(obj.xcur.thresvmap(:)),'%.2f') ' / ' num2str(obj.xcur.thresv0,'%.2f')]);
    end
    
    function showCutoff(obj,j,~)
      obj.loadframe(j);
      x = obj.xcur.x;
      m = false(size(obj.xcur.r2b));
      m(obj.xcur.neighInfo) = true;
      x = vmlColorify(x,m,2,64);
      x = vmlColorify(x,~obj.mfi.sm,1:3,-255);
      image(obj.mfi.yy,obj.mfi.xx,x); axis off;
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5 ]);
    end
    
    function showR2BCutoff(obj,j,~)
      obj.loadframe(j);
      x = repmat(uint8(obj.xcur.r2b*192),[1 1 3]);
      m = false(size(obj.xcur.r2b));
      m(obj.neighInfo(abs(diff(obj.xcur.r2b(obj.neighInfo),...
        [],2))>=obj.conf.thres.r2b_cutoff,:)) = true;
      x = vmlColorify(x,logical(m),2,64);
      x = vmlColorify(x,~obj.mfi.sm,1:3,-255);
      image(obj.mfi.yy,obj.mfi.xx,x); axis off;
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5 ]);
    end
    
    function showLumCutoff(obj,j,~)
      obj.loadframe(j);
      x = repmat(uint8(obj.xcur.lum*192),[1 1 3]);
      m = false(size(obj.xcur.r2b));
      m(obj.neighInfo(abs(diff(obj.xcur.lum(obj.neighInfo),...
        [],2))>=obj.conf.thres.lum_cutoff,:)) = true;
      x = vmlColorify(x,logical(m),2,64);
      x = vmlColorify(x,~obj.mfi.sm,1:3,-255);
      image(obj.mfi.yy,obj.mfi.xx,x); axis off;
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5 ]);
    end
    
    function showCutoff3(obj,j,thresv)
      clf;
      if nargin<3, thresv = []; end
      subplot(1,3,1); obj.showCutoff(j,thresv);
      subplot(1,3,2); obj.showR2BCutoff(j,thresv);
      subplot(1,3,3); obj.showLumCutoff(j,thresv);
    end

    function showCloudSeg(obj,j)
      obj.thres.v(:,j)=nan;
      obj.loadframe(0);
      clf;
      subplot(1,3,1);obj.showR2B(j);
      subplot(1,3,2);obj.showCutoff(j)
      subplot(1,3,3);obj.showThresVmap(j);
    end
    
    function showSobel(obj,j,~)
      x = obj.imread(j);
      x = vmlRed2Blue(x(:,:,1),x(:,:,3));
      x = vmlMedianBlur(uint8(x/max(x(:))*255),15);
      gb = 45;
      x = max(abs(vmlGaussianBlur(vmlSobelX(x),gb)),...
              abs(vmlGaussianBlur(vmlSobelY(x),gb)));
      x = vmlMinFilter(x,7);
      x(~obj.oi.sm) = NaN;
      assignin('base','X',x);
      imagesc(x); axis off; colorbar vert;
      %set(imagesc(x),'alphadata',~isnan(x));axis off;
    end
    
    function axlim = showPredIm(obj,j,tpred,zoom,linestyle)
      if nargin<4, zoom = 1; end
      if nargin<5, linestyle = 'g'; end
      mvec1 = obj.mvec_filtered(j);
      if any(isnan(mvec1)), obj.showframe(j); return; end
      
      yx = sunpos_plane(obj,j);
      mvnorm = sqrt(sum(mvec1.^2));
      mvec1 = mvec1/mvnorm;
      mvec1rot = [-mvec1(2);mvec1(1)];
      width = obj.conf.pred.width0+mvnorm*obj.conf.pred.width_v*tpred;
      tband = obj.conf.pred.tband0+obj.conf.pred.tband*tpred;
      
      tt = round(tpred-tband):round(tpred+tband);
      N = length(tt);
      tt1 = 0:round(tpred+tband);
      N1 = length(tt1);
      ww = (((0:N-1)/(N-1))*2-1)*width;
      
      if zoom<1, ymi = 1; xmi = 1; else
        axlim = do_lines([]); 
        xmi = axlim(1); ymi = axlim(3);
      end
      if zoom>=0
        x = obj.imread(j);
        if zoom>=1, x = x(axlim(3):axlim(4),axlim(1):axlim(2),:); end
        image(x); axis off;        
      end
      hold on;
      axlim = do_lines(linestyle);
      hold off;
      title([datestr(obj.ti(j),'HH:MM:SS') ' #' num2str(j) ' tpred=' num2str(tpred)]);

      function axlim = do_lines(linestyle)
        axlim = [inf -inf inf -inf];
        for i=-(zoom<2):obj.conf.pred.ngrid+(zoom<2)
          if zoom<2 && any(i==(obj.conf.pred.ngrid-1)/2+[0 1])
            yx1 = obj.camworld2im([repmat(yx+mvec1rot*width*...
              (2*i/obj.conf.pred.ngrid-1),1,N1)-mvec1*mvnorm*tt1;ones(1,N1)]);
          else
            yx1 = obj.camworld2im([repmat(yx+mvec1rot*width*...
              (2*i/obj.conf.pred.ngrid-1),1,N)-mvec1*mvnorm*tt;ones(1,N)]);
          end
          axlim([1 3]) = min(axlim([1 3]),fliplr(floor(min(yx1,[],2)')));
          axlim([2 4]) = max(axlim([2 4]),fliplr(ceil(max(yx1,[],2)')));
          if ~isempty(linestyle) && i>=0 && i<=obj.conf.pred.ngrid
            plot(yx1(2,:)-xmi+1,yx1(1,:)-ymi+1,linestyle); 
          end
          yx1 = obj.camworld2im([repmat(yx-mvec1*mvnorm*(tpred+...
            tband*(2*i/obj.conf.pred.ngrid-1)),1,N)+mvec1rot*ww;ones(1,N)]);
          axlim([1 3]) = min(axlim([1 3]),fliplr(floor(min(yx1,[],2)')));
          axlim([2 4]) = max(axlim([2 4]),fliplr(ceil(max(yx1,[],2)')));
          if ~isempty(linestyle) && i>=0 && i<=obj.conf.pred.ngrid
            plot(yx1(2,:)-xmi+1,yx1(1,:)-ymi+1,linestyle); 
          end
        end
      end
    end
    
    function axlim = showPredIm0(obj,j,tpred,linestyle)
      if nargin<4, linestyle = 'g'; end
      axlim = obj.showPredIm(j,tpred,0,linestyle);
    end
    
    function plotfeatures(obj,j,tpred)
      if nargin<3, tpred = []; end
      Z1 = obj.Z{obj.getZidx(tpred)};
      visualize_ngrids(Z1(:,j),obj.conf.pred.n_feature_maps,0,1);
      nstr = num2str(obj.conf.pred.nsubgrid);
      title(['red/blue ratio, mean/min/max over the median of a ' nstr 'x' nstr ' subgrid']);
    end
    
    function show(obj, tpred, method)
      if nargin<2, tpred = 0; end
      if nargin<3, method = 'showtraj'; end
      obj.newfig(1327,['PV power, tpred = ' num2str(tpred)]);
      subplot(3,2,5:6);
      obj.plotpow(tpred);
      title(obj.strtitle,'interpreter','none');
      j = find(all(~isnan(obj.mvec),1),1);
      if isempty(j), j=1; end
      show_upd(j);
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@show_upd)
      function out = show_upd( j, event_obj)
        if nargin>=2
          [~,j] = min(abs(event_obj.Position(1)-obj.ti));
        end
        figure(1327);
        subplot(3,2,1:4); 
        obj.(method)(j,tpred);
        drawnow;
        if nargin<2, out=[]; 
        else out = datestr(event_obj.Position(1)+tpred/86400,'HH:MM:SS'); end
      end
    end
    
    function show3(obj, tpred, method1, method2, method3)
      method2default = 'showR2B';
      method3default = 'showLum';
      if nargin<2, tpred = 0; end
      if nargin<3, method1 = 'showPredIm0';
      elseif isnumeric(method1)
        switch method1
          case 0, method1 = 'showPredIm0';
          case 1 
            method1 = 'showPredIm0';
            method2default = 'showRed2Blue100';
            method3default = 'plotThresholdScores';
          otherwise
            error('unknown code');
        end
      end
      if nargin<4, method2 = method2default; end
      if nargin<5, method3 = method3default; end
      obj.newfig(1325,['3 plots, tpred = ' num2str(tpred)]);
      subplot(2,3,4:6);
      obj.plotpow(tpred);
      title(obj.strtitle,'interpreter','none');
      j = find(all(~isnan(obj.mvec),1),1);
      if isempty(j), j=1; end
      show3_upd(j);
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@show3_upd)
      function out = show3_upd( j, event_obj)
        if nargin>=2
          [~,j] = min(abs(event_obj.Position(1)-obj.ti));
        end
        figure(1325);
        hax(1) = subplot(2,3,1); 
        obj.(method1)(j,tpred);
        hax(2) = subplot(2,3,2); 
        obj.(method2)(j,tpred);
        hax(3) = subplot(2,3,3); 
        obj.(method3)(j,tpred);
        xlim = get(hax,'xlim');
        ylim = get(hax,'ylim');
        linkaxes(hax([2;find(all(diff(cat(1,xlim{:}))==0,2) & ...
                             all(diff(cat(1,ylim{:}))==0,2))*2-1]))
        drawnow;
        if nargin<2, out=[]; 
        else out = datestr(event_obj.Position(1)+tpred/86400,'HH:MM:SS'); end
      end
    end
    
    function showmvec(obj)
      if all(isnan(obj.mvec(:))), error('no motion vectors available'); end
      obj.newfig(1326,'Motion vectors');
      subplot(3,2,5:6);
      v = obj.mvec_filtered;
      plot(obj.ti,obj.mvec(1,:),'-b.',obj.ti,obj.mvec(2,:),'-r.',...
        obj.ti,v(1,:),'b',obj.ti,v(2,:),'r');
      datetickzoom;
      title(obj.strtitle,'interpreter','none');
      j = find(all(~isnan(obj.mvec),1),1);
      if isempty(j), j=1; end
      show_mvec_upd(j);
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@show_mvec_upd)
      function out = show_mvec_upd( j, event_obj)
        use_filt = 1;
        if nargin>=2
          [~,j] = min(abs(event_obj.Position(1)-obj.ti));
          if all(~isnan(obj.mvec(:,j)))
            [~,j2] = min(abs([obj.mvec(:,j);obj.mvec_filtered(j)]-event_obj.Position(2)));
            use_filt = floor((j2-1)/2);
          end
        end
        figure(1326);
        subplot(3,2,1:4); 
        obj.showtraj(j,0,use_filt);
        drawnow;
        if nargin<2, out=[]; 
        else out = datestr(event_obj.Position(1),'HH:MM:SS'); end
      end
    end
    
    function showpow3(obj)
      obj.newfig(1324,'PV power (3 plots)');
      subplot(2,3,4:6);
      obj.plotpow; hold on; 
      hmarks = plot(obj.ti(1:3),interp1(obj.tm,obj.P,obj.ti(1:3)),'ro'); 
      hold off;
      set(hmarks, 'XData', [], 'YData', []);
      title(obj.strtitle,'interpreter','none');
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@show3pow_upd)
      show3pow_upd(obj.ti(end));
      function out = show3pow_upd( dticlick, event_obj)
        if nargin>1, dticlick = event_obj.Position(1); end
        figure(1324);
        jj = zeros(1,3);
        for i=0:2
          [~,j] = min(abs(dticlick-i*obj.conf.show_dt/86400-obj.ti));
          subplot(2,3,3-i); 
          obj.showframe(j,1);
          jj(i+1) = j;
        end
        set(hmarks, 'XData', obj.ti(jj), 'YData', interp1(obj.tm,obj.P,obj.ti(jj)));
        drawnow;
        out = datestr(dticlick,'HH:MM:SS');
      end
    end
    
    function showfeatures(obj,tpred,method)
      if nargin<2, tpred = []; end
      if nargin<3, method = 'showPredIm'; end
      [~, tpred] = obj.getZidx(tpred);
      obj.newfig(1322,['Features, tpred = ' num2str(tpred)]);
      subplot(2,3,4:6);
      obj.plotpow(tpred);
      title(obj.strtitle,'interpreter','none');
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@showfeatures_upd)
      showfeatures_upd(obj.ti(2));
      function out = showfeatures_upd( dticlick, event_obj)
        if nargin>1, dticlick = event_obj.Position(1); end
        figure(1322);
        [~,j] = min(abs(dticlick-obj.ti));
        subplot(2,3,1:2); 
        obj.plotfeatures(j,tpred);
        subplot(2,3,3);
        if any(isnan(obj.mvec_filtered(j))), image(obj.imread(j)); axis off;
        else obj.(method)(j,tpred); end
        drawnow;
        out = datestr(dticlick+tpred/86400,'HH:MM:SS');
      end
    end
    
    function showfea0(obj,tpred)
      if nargin<2, tpred = []; end
      obj.showfeatures(tpred,'showPredIm0');
    end

    function showcomp(obj,t1,t2,PBP)
      if nargin<4, PBP=[]; end
      j = obj.getidx(t1,t2);
      tt=obj.ti(j(1)):1/86400:obj.ti(j(end));
      plot(tt,obj.getP(tt+120/86400),'b',obj.ti(j),obj.pred(j,PBP),'r');
      datetickzoom;
    end
        
% end of plotting functionality
    
%% image <-> camera world transformations
    
    function p = im2camworld(obj,yx,yx0)
      %image pixel coordinates to camera world
      %yx0 is the internal offset (due to sky_area)
      if nargin<3, yx0 = obj.oi.sky_area([1 3])-1; end
      yx0 = yx0+obj.conf.img_offset;
      p = cam2world(yx+repmat(yx0(:),1,size(yx,2)),obj.calib.intern);
      p(3,:) = -p(3,:);
    end
    
    function yx = camworld2im(obj,p,yx0)
      %camera world to image pixel coordinates
      %yx0 is the internal offset (due to sky_area)
      if nargin<3, yx0 = obj.oi.sky_area([1 3])-1; end
      yx0 = yx0+obj.conf.img_offset;
      p(3,:) = -p(3,:);
      yx = world2cam(p,obj.calib.intern)-repmat(yx0(:),1,size(p,2));
    end
    
% end of image <-> camera world transformations
    
%% position of the sun
    function calc_sun_positions(obj)
        pos = zeros(3,length(obj.ti));
        [center, zenith] = sunpos_midday(obj);
        use_corrected_pos = zeros(1,length(obj.ti));
        pos(1,:) = obj.ti;
        for j=1:length(obj.ti)
            [Q, sun_zenith] = obj.sunpos_realworld_v2(j);
            pos([2 3],j) = obj.camworld2im(obj.ext_calib.R'*Q);
            if (sun_zenith-zenith)*zenith > 710 % this formula is empricially calculated to only consider a band close to image borders for correction
                use_corrected_pos(j) = 1;
            end
        end
         corrected_pos = obj.adjust_sun_pos(pos(2:3,:));
         pos(2:3,:) = corrected_pos.*[use_corrected_pos;use_corrected_pos] + pos(2:3,:).*[1-use_corrected_pos;1-use_corrected_pos]; 
         obj.sun_pos = pos;
    end
    function inv_rot_pos = adjust_sun_pos(obj,pos)
        [center, min_zenith,~] = sunpos_midday(obj);
        sh = .033; % the shift parameter and rotation angle is emprically set, but can be relate to min_zenith of day too.
        theta = degtorad(130); % rotaion angle counterclockwise about the center
        Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        Rot_inv = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
        
        %center = fix((pos(:,1)+pos(:,end))/2);
        center_pos = repmat(center, 1, size(pos,2));
        % shift points in the plane so that the center of rotation is at the origin
        rot_pos = Rot*(pos - center_pos) + center_pos;
        [min_y, idx] = min(rot_pos(2,:));
        min_x = rot_pos(1,idx);
        rot_pos(2,:) = (rot_pos(2,:)-min_y).*1+(sh^2) + min_y;
        rot_pos(1,:) = (rot_pos(1,:)-min_x).*(1-sh) + min_x;
        inv_rot_pos = Rot_inv*(rot_pos - center_pos) + center_pos;
    end
    function inv_rot_pos = adjust_sun_pos_v2(obj,pos)
        % works with real world 3d sun position data
        [~, min_zenith, center] = sunpos_midday(obj);
        center = center(1:2);
        sh = .033; % the shift parameter and rotation angle is emprically set, but can be relate to min_zenith of day too.
        theta = degtorad(130); % rotaion angle counterclockwise about the center
        Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        Rot_inv = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
        
        %center = fix((pos(:,1)+pos(:,end))/2);
        center_pos = repmat(center, 1, size(pos,2));
        % shift points in the plane so that the center of rotation is at the origin
        rot_pos = Rot*(pos - center_pos) + center_pos;
        [min_y, idx] = min(rot_pos(2,:));
        min_x = rot_pos(1,idx);
        rot_pos(2,:) = (rot_pos(2,:)-min_y).*1+(sh^2) + min_y;
        rot_pos(1,:) = (rot_pos(1,:)-min_x).*(1-sh) + min_x;
        inv_rot_pos = Rot_inv*(rot_pos - center_pos) + center_pos;
    end
    
    function t= gen_time_struct(obj,j,t_array)
      t = [];
      [t.year,t.month,t.day,t.hour,t.min,t.sec]=datevec(obj.ti(j));
      if exist('t_array','var')
          t.hour = t_array(1);
          t.min = t_array(2);
          t.sec = t_array(3);       
      end
      t.UTC = obj.calib.model3D.UTC;
      % for Daylight saving
      if ~isdst(datetime(t.year,t.month,t.day,t.hour,t.min,t.sec,'TimeZone',obj.conf.timezone))
        t.hour = t.hour + 1;
      end  
    end
    function p = sunpos_realworld(obj,j)
      %position of the sun in the world coordinates
      t = obj.gen_time_struct(j);
      sun = sun_position(t, obj.calib.model3D.Location);
      [x,y,z] = sph2cart(sun.azimuth*pi/180,(90-sun.zenith)*pi/180,1);
      p = [x y z]';
    end
    function [p, sun_zenith] = sunpos_realworld_v2(obj,j)
      %position of the sun in the world coordinates
      t = obj.gen_time_struct(j);
      sun = sun_position(t, obj.calib.model3D.Location);
      [x,y,z] = sph2cart(sun.azimuth*pi/180,(90-sun.zenith)*pi/180,1);
      p = [x y z]';
      sun_zenith = sun.zenith;
    end
    
    function [p, sun_zenith] = sunpos_realworld_v3(obj,tm)
      %position of the sun in the world coordinates
      t = [];
      [t.year,t.month,t.day,t.hour,t.min,t.sec]=datevec(tm);
      t.UTC = obj.calib.model3D.UTC;
      % for Daylight saving
      if ~isdst(datetime(t.year,t.month,t.day,t.hour,t.min,t.sec,'TimeZone',obj.conf.timezone))
        t.hour = t.hour + 1;
      end  
      sun = sun_position(t, obj.calib.model3D.Location);
      [x,y,z] = sph2cart(sun.azimuth*pi/180,(90-sun.zenith)*pi/180,1);
      p = [x y z]';
      sun_zenith = sun.zenith;
    end
    
    function yx = sunpos_im(obj,j)
      %position of the sun in image pixels
      yx = obj.camworld2im(obj.calib.Rext'*obj.sunpos_realworld(j));
    end
    
    function p = sunpos_plane(obj,j)
      %position of the sun in the unit plane
      p = obj.calib.Rext'*obj.sunpos_realworld(j);
      p = p(1:2)/p(3);
    end
    
    function P = plant_projection_on_image(obj,j,h)
        % project plant shape onto sky image
        m = obj.sunpos_realworld(j); % line slope for projection
        m(1:2)=obj.adjust_sun_pos_v2(m(1:2));
        if size(obj.conf.plant_coords,2)<1
            return;
        end
        P = obj.camworld2im(obj.calib.Rext'*(bsxfun(@plus,...
          obj.conf.plant_coords,m*h/m(3))));
        P(:,end+1) = P(:,1);
        P = obj.adjust_sun_pos(P);
    end


% end position of the sun
    
%% sky masks
    
    function sm = skymask_wo_sun(obj,j)
      %sky mask without sun, frame j (taking out a circle around 
      %the sun position)
      yx = sunpos_im(obj,j);
      sm = obj.oi.sm;
      sm((obj.oi.XX-yx(2)).^2+(obj.oi.YY-yx(1)).^2<=obj.conf.rsun_px^2) = false;
    end
    
    function SM = skymask_oflow(obj,j,mvec1)
      %compute sky mask for the part of the sky where the clouds
      %come from
      mvec1 = mvec1/sqrt(sum(mvec1.^2));
      sm = obj.skymask_wo_sun(j);
      yx = sunpos_plane(obj,j);
      mvec1rot = [-mvec1(2);mvec1(1)];
      llim = yx'*mvec1rot-obj.conf.oflow.mask_width; 
      hlim = yx'*mvec1rot+obj.conf.oflow.mask_width;
      if obj.conf.oflow.mask_cover_center
        if llim<0 && hlim<0, hlim=0; elseif llim>0 && hlim>0, llim=0; end
      end
      ppxproj = obj.oi.ppx*mvec1rot;
      sm(ppxproj<llim | ppxproj>hlim) = false;
      sm_upwind = sm;
      sm_upwind(obj.oi.ppx*mvec1-yx'*mvec1>0) = false;
      SM = cat(3,sm,sm_upwind);
    end
      
    
% end sky masks
    
%% sun detection / external calibration (Kabsch algorithm)
    
    function x = im4sundetect(obj,j)
      %compute image for automatic detection of the sun
      x = sum(obj.imread(j),3)/(3*255);
      x(~(obj.oi.sm | obj.oi.sm_border)) = 0.5;
      x = vmlSobelAbs(x);
      x(~obj.oi.sm) = 0;
      x=vmlGaussianBlur(x,obj.conf.gaussblur4sundect_px);
    end
    
    function yx = detect_saturated_sun(obj,j)
      %automatic detection of the (saturated) sun
      x = obj.im4sundetect(j);
      [~,jmax] = max(x(:));
      [p,q] = ind2sub(obj.sz,jmax);
      yx = [p;q];
      %figure(55);imagesc(x);hold on;plot(q,p,'gx');hold off
    end
    
    function sat_sun_pos=detect_saturated_sun_v2(obj,j)
      sat_sun_pos=[];
      sun_pos = obj.sun_pos([2 3],j);
      sat_sun_pos_old = obj.detect_saturated_sun(j);
      if norm(sat_sun_pos_old-sun_pos)<30
        sun_pos = sat_sun_pos_old;
      end
      
      eccthresh = 0.6;
      areaThresh = 600;
      
      cx=floor(sun_pos(1)); cy=floor(sun_pos(2)); R=30;
      im = obj.imread(j);
      
      sun_p = im(max(cx-R,0):min(cx+R,size(im,1)),max(cy-R,0):min(cy+R,size(im,2)),:); 
      
      xg=rgb2gray(sun_p);
      xg=vmlGaussianBlur(xg,obj.conf.gaussblur4sundect_px);
      gthres=0.4*mean(mean(xg));
      xg2=xg;
      xg(xg2>gthres)=0;
      xg(xg2<=gthres)=255;

        %Boundary Label the Filtered Image
        [L, num]=bwlabel(xg);
        STATS=regionprops(L,'all');
        removed=0;
        %Remove the noisy regions 
        for j=1:num
            dd=STATS(j).Area;
            if (dd < 5 || (dd > areaThresh && STATS(j).Eccentricity > eccthresh))
                L(L==j)=0;
                removed = removed + 1; 
            end
        end
        if removed>0
            [L2, num]=bwlabel(L);
            STATS=regionprops(L2,'all');
        end
        max_area=0;
        for j=1:num
            if STATS(j).Area > max_area
                sat_sun_pos = STATS(j).Centroid + [cy-R,cx-R];
                max_area=STATS(j).Area;
            end
        end
        
    end
    
    function R = kabsch(obj,R,frame_nos)
      %Kabsch algorithm for compution rotation matrix such that
      %detected sun position and acutal sun position coincide
      if nargin<2, R = []; end
      if nargin<3 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      pp = zeros(2,length(frame_nos));
      PP = zeros(3,length(frame_nos)); QQ = PP;
      for j = 1:length(frame_nos)
        pp(:,j) = obj.detect_saturated_sun(frame_nos(j))';
        PP(:,j) = obj.im2camworld(pp(:,j));
        QQ(:,j) = obj.sunpos_realworld(frame_nos(j));
      end
      if isempty(R)
        [U,~,V] = svd(PP*QQ');
        R = V*diag([1 1 sign(det(V*U'))])*U';
      end
      obj.calib.kabsch.frame_nos = frame_nos;
      obj.calib.kabsch.dt = obj.ti(frame_nos);
      obj.calib.kabsch.world_detect = PP;
      obj.calib.kabsch.realword = QQ;
      obj.calib.kabsch.im_detect = pp;
      obj.calib.kabsch.im_model = obj.camworld2im(R'*QQ);
      obj.calib.kabsch.R = R;
      if obj.printlevel>0
        disp(['#frames = ' num2str(length(frame_nos)) ', pixel RMSE = ' ...
          num2str(sqrt(mean(sum((pp-obj.calib.kabsch.im_model).^2,1))),'%.1f') ...
          ', RMSE on 1 m unit sphere = ' ...
          num2str(1000*sqrt(mean(sum((R*PP-QQ).^2,1))),'%.1f') ' mm']);
      end
    end
    
% end sun detection

%% motion vectors / optical flow

    function remove_zero_mvec(obj, tol, frame_nos)
      if nargin<3 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      r = false(1,length(obj.ti));
      r(frame_nos) = (sqrt(sum(obj.mvec(:,frame_nos).^2,1))<tol);
      obj.mvec(:,frame_nos(r)) = NaN;
    end
    
    function v = mvec_filtered(obj, frame_nos)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end 
      v = zeros(2,numel(frame_nos)) + NaN;
      for i = 1:numel(frame_nos)
        j = frame_nos(i);
        if j<2 || j>length(obj.ti), continue; end
        if j-1>=obj.conf.nfilter_short_mvec
          n = obj.conf.nfilter_short_mvec;
          vv = obj.mvec(:,j-n+1:j);
          if all(~isnan(vv(:)))
            bb = zeros(2);
            for ii=1:2
              bb(ii,:) = [ones(1,n);(1:n)-obj.conf.evfilter_short_mvec]'\vv(ii,:)';
            end
            stdev = max(sqrt(mean((repmat(bb(:,1),1,n)+repmat(bb(:,2),1,n).*...
              repmat((1:n)-obj.conf.evfilter_short_mvec,2,1)-vv).^2,2)))/max(norm(bb(:,1)),1e-8);
            if stdev<=obj.conf.filter_short_maxstd
              v(:,i) = bb(:,1);
              continue;
            end
          end
        end
        n = min(obj.conf.nfilter_mvec,j-1);
        vv = obj.mvec(:,j-n+1:j);
        k = find(all(isnan(vv),1),1,'last');
        if isempty(k), k=0; end
        if k<n 
          v(:,i) = mean(vv(:,k+1:end),2);
        else
          k = find(all(~isnan(obj.mvec(:,1:j)),1),1,'last');
          if isempty(k), v(:,i) = NaN; else v(:,i) = obj.mvec(:,k); end
        end
      end
    end
    
    function mvec1 = optimize_optical_flow(obj, j, SM, mvec_start, diagnose)
      if nargin<5, diagnose=0; end
      x0 = obj.imread(j-1);
      x0 = double(vmlMedianBlur(x0(:,:,1),obj.conf.oflow.median_ksize));
      x1 = obj.imread(j);
      x1 = double(vmlMedianBlur(x1(:,:,1),obj.conf.oflow.median_ksize));
      
      sm = SM(:,:,1);
      if size(SM,3)==2, sm_upwind = SM(:,:,2);
      else sm_upwind = false(size(sm)); end
      jsm = find(sm);
      jsm1 = jsm;

      dx = [0 0];
      dt = (obj.ti(j)-obj.ti(j-1))*86400;
      dx1 = dx; jsm11 = jsm;
      evaluate_fit4oflow_count = diagnose-1;
      fit = evaluate_fit4oflow;
      
      for i = 1:size(mvec_start,2)
        dx1 = round((-mvec_start(:,i)'*dt)/obj.dplane)*obj.dplane;  
        jsm11=closest_in_plane(obj.oi.ppx,jsm,jsm1,dx1,size(sm),0);
        fit1 = evaluate_fit4oflow;
        if fit1<fit, dx = dx1; jsm1 = jsm11; fit = fit1; end
      end
      
      evaluate_fit4oflow_count = [evaluate_fit4oflow_count -1];
      scale = 25;
      iprev = [];
      while 1
        j1best = [];
        for i=0:3
          if ~any(i==iprev)
            dx1 = dx + scale*obj.dplane*(-1)^rem(i,2)*([0 1]+floor(i/2)*[1 -1]);
            jsm11=closest_in_plane(obj.oi.ppx,jsm,jsm1,dx1,size(sm),0+(scale==1));
            fit1 = evaluate_fit4oflow;
            if fit1<fit
              dx1best = dx1; fit = fit1; j1best = jsm11; iprev1 = bitxor(i,1);
            end
          end
        end
        if ~isempty(j1best)
          dx = dx1best; jsm1 = j1best; iprev = iprev1;
          if obj.printlevel+2*any(evaluate_fit4oflow_count>=0)>=3
            disp(['*** ' sprintf('%3i',round(dx/obj.dplane)) ' ' sprintf('%12i',fit)]);
          end
        else
          if scale==1, break; end
          iprev = [];
          scale = round(scale/5);
        end
      end
      if evaluate_fit4oflow_count(1)>=0
        evaluate_fit4oflow_count = evaluate_fit4oflow_count(1);
        dx1 = dx; jsm11 = jsm1; evaluate_fit4oflow;
      end
      mvec1 = -dx(:)/dt;
      
      function fit1 = evaluate_fit4oflow
        sm11 = false(size(sm)); sm11(jsm11) = true;
        diffx=zeros(size(sm)); 
        diffx(jsm) = abs(x0(jsm11)-x1(jsm));
        diffx(~(sm & sm11))=nan;
        fit1 = sum((1+obj.conf.oflow.w_upwind*sm_upwind(~isnan(diffx))).*...
               max(0,diffx(~isnan(diffx))-obj.conf.oflow.ignore_diff).^2);
        if obj.printlevel+2*any(evaluate_fit4oflow_count>=0)>=3
          if all(evaluate_fit4oflow_count>=0), spc = '+++ ';
          else spc = '    '; end
          disp([spc sprintf('%3i',round(dx1/obj.dplane)) ' ' sprintf('%12i',fit1)]);
        end
        if obj.printlevel>=1 && all(evaluate_fit4oflow_count>=0)
          figure(200+evaluate_fit4oflow_count);
          set(imagesc(diffx),'alphadata',~isnan(diffx));axis off;
          colorbar vert;
          assignin('base',['XX' num2str(evaluate_fit4oflow_count)],diffx);
          if evaluate_fit4oflow_count==0, assignin('base','SM_upwind',sm_upwind); end
          evaluate_fit4oflow_count = evaluate_fit4oflow_count+1;
        end
      end
    end
    
    function optical_flow(obj, frame_nos, do_replace)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 2:length(obj.ti); 
      end
      if nargin<3, do_replace = 0; end
      for j = frame_nos(:)'
        if j<2 || j>length(obj.ti) || (~do_replace && all(~isnan(obj.mvec(:,j))))
          continue; 
        end
        if obj.printlevel>=1
          disp(['computing motion vector for frame #' num2str(j)]);
        end
        
        if j<=3, mvec_start = obj.mvec(:,j-1); 
        else mvec_start = [obj.mvec(:,j-1) obj.mvec_filtered(j-1)]; end
        mvec_start(:,any(isnan(mvec_start),1)) = [];
        
        if isempty(mvec_start) || ~obj.conf.oflow.masked
          sm = obj.skymask_wo_sun(j);
          mvec1 = optimize_optical_flow(obj, j, sm, mvec_start, length(frame_nos)==1 && do_replace);
          mvec_start = mvec1;
        end
        
        if obj.conf.oflow.masked
          SM = obj.skymask_oflow(j, mvec_start(:,end));
          mvec1 = optimize_optical_flow(obj, j, SM, mvec_start, length(frame_nos)==1 && do_replace);
        end
        obj.mvec(:,j) = mvec1;
      end
    end
  
% end motion vectors

%% cloud thresholding

    function opt_cur_thresv(obj)
      obj.xcur.thresv0 = mean(obj.xcur.neighy);
      method = 'simple'; %simple, wavg, qp
      if strcmp(method,'simple') || (strcmp(method,'wavg') && numel(obj.thres.wavg.yy)==1)
        obj.xcur.thresv = obj.xcur.thresv0;
        obj.xcur.thresvmap = zeros(size(obj.mfi.sm))+obj.xcur.thresv;
      elseif strcmp(method,'wavg')
        %weighted averaging with Gaussian Kernel
        K=gausskern(obj.thres.wavg.sigma,obj.xcur.neighp',...
          [obj.thres.wavg.xx(:)';obj.thres.wavg.yy(:)']);
        K = K + obj.conf.thres.wavg.offset;
        obj.xcur.thresv = reshape(1./sum(K).*sum(bsxfun(@times,K,...
          obj.xcur.neighy)),size(obj.thres.wavg.yy));
        obj.xcur.thresvmap = interp2(obj.thres.wavg.xx,obj.thres.wavg.yy,...
          obj.xcur.thresv,obj.thres.XX,obj.thres.YY);
      else
        %optimized threshold map taking into account same class and
        %opposite class neighbor pairs
        means = nan(obj.conf.thres.ngrid);
        nmeans = zeros(obj.conf.thres.ngrid);
        plateaus = nan(obj.conf.thres.ngrid^2,4);
        k = 1;
        for i = 1:obj.conf.thres.ngrid
          for j = 1:obj.conf.thres.ngrid
            jj = (obj.xcur.neighp(:,1)>=obj.thres.xxg(i) & ...
              obj.xcur.neighp(:,1)<obj.thres.xxg(i+1) & ...
              obj.xcur.neighp(:,2)>=obj.thres.yyg(j) & ...
              obj.xcur.neighp(:,2)<obj.thres.yyg(j+1));
            nmeans(k) = nnz(jj);
            means(k) = mean(obj.xcur.neighy(jj));
            jj = (obj.xcur.nneighp(:,1)>=obj.thres.xxg(i) & ...
              obj.xcur.nneighp(:,1)<obj.thres.xxg(i+1) & ...
              obj.xcur.nneighp(:,2)>=obj.thres.yyg(j) & ...
              obj.xcur.nneighp(:,2)<obj.thres.yyg(j+1));
            if nnz(jj)>=200
              v0 = means(k);
              if isnan(v0), v0 = obj.xcur.thresv0; end
              plateaus(k,:) = find_plateau(obj.xcur.nneighy(jj),v0,50,5,5);
            end
            k=k+1;
          end
        end
        obj.xcur.thresv=SmoothMapQP(means,nmeans,size(obj.neighInfo,1)/1000,0,1,0.1,plateaus);
        VVg = obj.thres.XXg;
        VVg(2:end-1,2:end-1) = obj.xcur.thresv;
        VVg(1,:) = VVg(2,:); VVg(end,:) = VVg(end-1,:);
        VVg(:,1) = VVg(:,2); VVg(:,end) = VVg(:,end-1);
        obj.xcur.thresvmap = interp2(obj.thres.XXg,obj.thres.YYg,VVg,...
            obj.thres.XX,obj.thres.YY);
      end
      obj.xcur.r2bthres(~obj.mfi.sm) = NaN;
    end

% end cloud thresholding

%% features / image parts for prediction

    function tp = t_predictable(obj,frame_nos)
      if nargin<3 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      tp = zeros(1,length(frame_nos))+NaN;
      for i = 1:length(frame_nos)
        j = frame_nos(i);
        mvec1 = obj.mvec_filtered(j);
        if all(~isnan(mvec1))
          yx = obj.sunpos_im(j);
          p0 = obj.sunpos_plane(j);
          if t_predictable_f(900)<=0, tp(i) = inf;
          else tp(i) = fzero(@t_predictable_f,[0 900]); end
        end
      end
      function y = t_predictable_f(t)
        y = sqrt(sum((obj.camworld2im([p0-mvec1*t;1])-yx).^2))-obj.conf.rsun_px;
      end
    end
    
    function [ZZ, npoints1] = calc_features(obj, tpred, j, N, Nsub)
      if nargin<4, N = obj.conf.pred.ngrid; end
      if nargin<5
        if N == obj.conf.pred.ngrid, Nsub = obj.conf.pred.nsubgrid;
        else Nsub = 1; end
      end
      mvec1 = obj.mvec_filtered(j);
      if any(isnan(mvec1)), error('no motion vector available'); end
      
      yx = sunpos_plane(obj,j);
      mvnorm = sqrt(sum(mvec1.^2));
      mvec1 = mvec1/mvnorm;
      mvec1rot = [-mvec1(2);mvec1(1)];
      width = obj.conf.pred.width0+mvnorm*obj.conf.pred.width_v*tpred;
      tband = obj.conf.pred.tband0+obj.conf.pred.tband*tpred;
      
      x = double(reshape(obj.imread(j),numel(obj.oi.sm),3));
      x = x(obj.oi.sm,[1 3])/255;
      ppx = obj.oi.ppx(obj.oi.sm,:);
      jj = abs(ppx*mvec1rot-yx'*mvec1rot)>width | ...
        abs(ppx*mvec1-yx'*mvec1+tpred*mvnorm)>tband*mvnorm;
      ppx(jj,:) = []; x(jj,:) = [];
      
      k = 1;
      npoints1 = [0;inf];
      subgrid = zeros(Nsub);
      ZZ = cell(1,obj.conf.pred.n_feature_maps);
      for i=1:obj.conf.pred.n_feature_maps, ZZ{i} = zeros(N); end
      for ix=1:N
        for iy=1:N
          for isubx=1:Nsub
            for isuby=1:Nsub
              jj = abs(ppx*mvec1rot-yx'*mvec1rot+width*...
                ((2*iy-1)/N-1+((2*isuby-1)/Nsub-1)/N))<=width/(N*Nsub) & ...
                abs(ppx*mvec1-yx'*mvec1+(tpred+tband*...
                ((2*ix-1)/N-1+((2*isubx-1)/Nsub-1)/N))*mvnorm)<=tband*mvnorm/(N*Nsub);
              subgrid(isubx,isuby) = median(vmlRed2Blue(x(jj,1),x(jj,2)));
              npoints1(1) = npoints1(1)+nnz(jj);
              npoints1(2) = min(npoints1(2),nnz(jj));
            end
          end
          ZZ{1}(k) = mean(subgrid(:));
          ZZ{2}(k) = min(subgrid(:));
          ZZ{3}(k) = max(subgrid(:));
          k = k+1;
        end
      end
      npoints1(1) = npoints1(1)/(N*Nsub)^2;
      for i=1:obj.conf.pred.n_feature_maps, ZZ{i} = flipud(ZZ{i}); end
    end
    
    function features(obj, tpred, frame_nos, do_replace)
      if nargin<2, error('please specify tpred'); end
      if nargin<3 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      if nargin<4, do_replace = 0; end
      N = obj.conf.pred.ngrid;
      jZ = obj.getZidx(tpred, 1);
      if isempty(jZ)
        Z1 = nan(2*N^2,length(obj.ti));
        npoints = nan(2,length(obj.ti));
        jZ = length(obj.Z_tpred)+1;
      elseif length(jZ)==1
        Z1 = obj.Z{jZ};
        npoints = obj.Z_npoints{jZ};
      else
        error('duplicate entry found in Z_tpred');
      end
      for j = frame_nos(:)'
        if any(isnan(obj.mvec_filtered(j))) || (all(~isnan(Z1(:,j))) && ~do_replace)
          continue;
        end
        if obj.printlevel>=1
          disp(['computing feature vector for frame #' num2str(j)]);
        end
        [ZZ, npoints1] = obj.calc_features(tpred, j);
        Z1(:,j) = reshape(cat(3,ZZ{:}),size(Z1,1),1);
        npoints(:,j) = npoints1;
        obj.Z{jZ} = Z1;
        obj.Z_npoints{jZ} = npoints;
        obj.Z_tpred(jZ) = tpred;
      end
    end
        
% end features / image parts for prediction

%% SVM
    function yy = persist_cloudiness(obj,tpred,frame_no_from,frame_no_to)
      if nargin<3, frame_no_from = 1; end
      if nargin<4, frame_no_to = length(obj.ti); end
      jj = frame_no_from:frame_no_to;
      yy = obj.cloudiness_label(obj.ti(jj));
      Yte = obj.cloudiness_label(obj.ti(jj) + tpred/86400);
      if obj.printlevel>=1
        fprintf('persistance classification error, frame %i to %i\n',frame_no_from,frame_no_to);
        disp_svm_err(Yte,yy);
      end
    end

    function [w,b] = trainsvm(obj,tpred,frame_nos)
      if nargin<3, frame_nos = 1; end
      gradb = obj.conf.svm.gradb_train;
      
      Z1 = obj.Z{obj.getZidx(tpred)};
      jtr = find(all(~isnan(Z1),1));
      if isscalar(frame_nos)
        jtr = jtr(1:round(length(jtr)*frame_nos));
      else
        jtr = intersect(jtr,frame_nos);
      end
      
      Z1 = Z1(:,jtr);
      Y = obj.cloudiness_label(obj.ti(jtr) + tpred/86400);
      j = find(Y~=0);
      d = size(Z1,1);
      N1 = sum(abs(Y));
      if gradb(1)
        nx = d+2*N1;
        nc = N1+2*(N1-1);
      else
        nx = d+1+N1;
        nc = N1;
      end

      Q = sparse(1:d,1:d,1,nx,nx);
      k = [zeros(1,nx-N1) zeros(1,N1)+obj.conf.svm.C];
      if gradb
        A1 = sparse([1:N1-1 1:N1-1],d+[1:N1-1 2:N1],[ones(1,N1-1) -ones(1,N1-1)],N1-1,nx);
        A = [-sparse(repmat(Y(j),d,1).*Z1(:,j))' sparse(1:N1,1:N1,-Y(j))...
          -speye(N1);A1;-A1];
        b = [-ones(N1,1);repmat(diff(j(:))*gradb,2,1)];
      else
        A = [-sparse([repmat(Y(j),d,1).*Z1(:,j);Y(j)])' -speye(N1)];
        b = -ones(nc,1);
      end
      lb = [-inf(1,nx-N1) zeros(1,N1)];
      ub = inf(1,nx);

      opts = [];
      opts.extractbounds = 0;
      opts.qfactor = 1;
      tic;
      [status,x] = ipoptqp(Q,k,A,b,[],[],lb,ub,opts);
      if obj.printlevel>=1
        fprintf('SVM trained for frames %i to %i, elapsed=%.1f s, status=%i\n',...
          min(jtr),max(jtr),toc,status);
      end

      w = x(1:d)';
      if gradb
        b = interp1(j,x(d+1:d+N1),1:length(jtr));
        j1 = find(~isnan(b),1); b(1:j1-1) = b(j1);
        j1 = find(~isnan(b),1,'last'); b(j1+1:end) = b(j1);
      else
        b = x(d+1);
      end

      if obj.printlevel>=1
        disp('training error');
        disp_svm_err(Y,saturate(w*Z1+b));
      end
    end
    
    function [yy,bb] = predsvm(obj,tpred,w,b,frame_no_from,frame_no_to)
      if nargin<5, frame_no_from = 1; end
      if nargin<6, frame_no_to = length(obj.ti); end
      jj = frame_no_from:frame_no_to;
      yy = nan(1,length(jj));
      bb = zeros(1,length(jj))+b;
      Z1 = obj.Z{obj.getZidx(tpred)};
      gradb = obj.conf.svm.gradb_pred;
      Yte = obj.cloudiness_label(obj.ti(jj) + tpred/86400);
      for i = 1:length(jj)
        j = jj(i);
        if any(isnan(Z1(:,j))), continue; end
        yy(i) = saturate(w*Z1(:,j)+bb(i));        
        bgoal =  bb(i) + abs(Yte(i))*(Yte(i)-yy(i));
        bb(obj.ti(jj)>=obj.ti(j) + tpred/86400) = ...
          max(bb(i)-gradb,min(bb(i)+gradb,bgoal));
      end
      if obj.printlevel>=1
        fprintf('SVM error, frame %i to %i\n',frame_no_from,frame_no_to);
        disp_svm_err(Yte,yy);
      end
    end

% end SVM
  end 
end
