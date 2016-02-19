classdef vmlSeq < handle
  properties
    data;           %image file names and measured data
    uidata;         %user interface data
    conf;           %configuration struct
    oi;             %original image info
    mfi;            %median filtered image info
    calib;          %camera calib
    calc;           %calculation results
    thres;          %data for cloud thresholding
    curseg;         %current segmentation result
    curmot;         %current motion estimation result
    curscale;       %current clear/cloudy to power / irr scaling
    curpred;        %current prediction
  end
  
  methods
  
%% basic functions: constructor, imread, getP, ...

    function obj = vmlSeq(plantid,folder,hours)
      if nargin<3, hours = []; end
      
      switch lower(plantid)
        case 'cavriglia'
          conf = vmlconfig_cavriglia;
        otherwise
          error(['unknown plantid ' plantid]);
      end
      obj.data.plantid = plantid;
      
      %each vmlSeq object gets an id, in order to identify it back
      %from Matlab figures
      obj.uidata.id = floor(rand*1e12);
      obj.uidata.printlevel = 1; %0=silent, 1=summarizing, 2=info, 3=verbose
      obj.uidata.varname = [];
      obj.uidata.open_figs = [];
      
      %store configuration and list image files
      obj.conf = conf;
      [obj.data.files, obj.data.folder, obj.data.ti] = vmlListImgFiles(folder, hours, conf);
      
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
      obj.calib.UTC = obj.calib.model3D.UTC+isdst(datetime(obj.data.ti(1),'ConvertFrom','datenum','TimeZone',conf.TimeZone));
            
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
      obj.oi.sz = size(obj.oi.sm);
      [obj.oi.YY,obj.oi.XX] = ndgrid(1:obj.oi.sz(1),1:obj.oi.sz(2));
      
      %compute border of the sky mask
      obj.oi.sm_border = vmlExtend1Px(obj.oi.sm) & ~obj.oi.sm;
      
      %reduce pixels on unit sphere to final bounding box and map
      %them to unit plane
      obj.oi.spx = obj.oi.spx(obj.oi.sky_area(1):obj.oi.sky_area(2),obj.oi.sky_area(3):obj.oi.sky_area(4),:);
      obj.oi.spx = reshape(obj.oi.spx,numel(obj.oi.sm),3);
      obj.oi.ppx = obj.oi.spx(:,1:2)./repmat(obj.oi.spx(:,3),[1 2]);
      obj.oi.ppx(~(obj.oi.sm | obj.oi.sm_border),:) = inf;
      
      %compute minimum distances of neighboring pixels on unit plane
      plane_px1 = reshape(obj.oi.ppx,[obj.oi.sz 2]);
      d1 = plane_px1(2:end,:,1)-plane_px1(1:end-1,:,1);
      d2 = plane_px1(:,2:end,2)-plane_px1(:,1:end-1,2);
      obj.oi.dplane = min(min(abs(d1(~isnan(d1)))),min(abs(d2(~isnan(d2)))));
      assert(obj.oi.dplane >= 1/max(height,width));
      
      %load irradiation and power data
      if exist([conf.tmpfolder folder '\data.mat'],'file')
        load([conf.tmpfolder folder '\data.mat']);
      else
        switch lower(plantid)
          case 'cavriglia'
            [data_Power, data_Power_ncols, data_Irr, data_Temp] = vml_loaddata_cavriglia(folder,conf);
        end
        save([conf.tmpfolder folder '\data.mat'],'data_Power', 'data_Power_ncols', 'data_Irr', 'data_Temp');
      end
      
      %crop irradiation and power data to the relevant time window
      obj.data.P = []; obj.data.Irr = []; obj.data.Temp = [];
      crop_fun = @(x)x(max(1,find(x(:,1)>=min(obj.data.ti),1)-1):min(size(x,1),...
        find(x(:,1)<=max(obj.data.ti)+conf.pad_minutes*60/86400,1,'last')+1),:);
      if ~isempty(data_Power), obj.data.P = crop_fun(data_Power); end
      if ~isempty(data_Irr), obj.data.Irr = crop_fun(data_Irr); end
      if ~isempty(data_Temp), obj.data.Temp = crop_fun(data_Temp); end
      obj.data.P_ncols = data_Power_ncols;
      obj.data.tmin = min(obj.data.ti);
      obj.data.tmax = max(obj.data.ti);
      
      if size(obj.data.P,1)>0
          %remove outliers in power
          count = 0;
          while count<3
            count = count+1;
            absdiff=abs(median([obj.data.P(1:end-2,2) obj.data.P(2:end-1,2) obj.data.P(3:end,2)],2)-obj.data.P(2:end-1,2));
            remove = [false absdiff'>mean(absdiff)+6*std(absdiff) false];
            if ~any(remove), break; end
            obj.data.P(remove,:) = [];
          end
      end
      
      %assign title and day in date time format
      obj.data.day = folder;
      obj.data.dt_day = datenum(obj.data.day,'yyyy_mm_dd');
      obj.data.day_of_year = 1+obj.data.dt_day-datenum(obj.data.day(1:4),'yyyy');
         
      %initialize data for motion vectors and power prediction
      obj.calc.mvec = nan(2,length(obj.data.ti));
      obj.calc.thres_v = nan(obj.conf.seg.ngrid^2,length(obj.data.ti));

      %do a median downscale of the image and compute the sky mask...
      [obj.mfi.sm,obj.mfi.xx,obj.mfi.yy,obj.mfi.ksize] = ...
        vmlMedianDownscale(obj.oi.sm,obj.conf.mfi.sz);
      obj.mfi.sm = logical(obj.mfi.sm);
      obj.mfi.sz = size(obj.mfi.sm);
      
      %... and the unit plane pixels for the median downscaled images 
      [obj.mfi.YY,obj.mfi.XX] = ndgrid(obj.mfi.yy,obj.mfi.xx);
      spx = obj.im2camworld([obj.mfi.YY(:)';obj.mfi.XX(:)'])';
      obj.mfi.ppx = spx(:,1:2)./repmat(spx(:,3),[1 2]);
      obj.mfi.ppx(~obj.mfi.sm,:) = nan;
      obj.mfi.ppx = vmlRemoveIsolatedPPX(obj.mfi);
      
      %compute a regular grid on the plane and the indexes of the closest
      %points in obj.mfi.ppx
      plane_px1 = reshape(obj.mfi.ppx,[obj.mfi.sz 2]);
      d1 = plane_px1(2:end,:,1)-plane_px1(1:end-1,:,1);
      d2 = plane_px1(:,2:end,2)-plane_px1(:,1:end-1,2);
      obj.mfi.dplane = min(min(abs(d1(~isnan(d1)))),min(abs(d2(~isnan(d2)))));
      wplane = max(obj.mfi.ppx)-min(obj.mfi.ppx);
      obj.mfi.prg_NN = round(wplane/obj.mfi.dplane)+2;
      obj.mfi.prg_yx0 = min(obj.mfi.ppx)+(wplane-(obj.mfi.prg_NN-1)*obj.mfi.dplane)/2;
      obj.mfi.prg_yy = obj.mfi.prg_yx0(1)+(0:obj.mfi.prg_NN(1)-1)*obj.mfi.dplane;
      obj.mfi.prg_xx = obj.mfi.prg_yx0(2)+(0:obj.mfi.prg_NN(2)-1)*obj.mfi.dplane;
      ppx = obj.mfi.ppx; ppx(isnan(ppx(:,1)),:) = 1e12;
      obj.mfi.prg_j = closest_in_plane(obj.mfi.prg_yy,obj.mfi.prg_xx,ppx,obj.mfi.sz);
      %check if every closest grid point brings back the point
      j = find(~isnan(obj.mfi.ppx(:,1)));
      jj = round(bsxfun(@minus,obj.mfi.ppx(j,:),obj.mfi.prg_yx0)/obj.mfi.dplane)+1;
      j2=obj.mfi.prg_j(jj(:,1)+obj.mfi.prg_NN(1)*(jj(:,2)-1));
      %accept presently one failure here
      assert(sum(j~=j2)<=1,'error in closest neighbor computation');
      
      %tiles for the local thresholding
      rr = ((1:obj.conf.seg.ngrid)-.5)/obj.conf.seg.ngrid;
      obj.thres.xxg = [1 round(1+(obj.mfi.sz(2)-1)*rr) obj.mfi.sz(2)]; 
      obj.thres.yyg = [1 round(1+(obj.mfi.sz(1)-1)*rr) obj.mfi.sz(1)]; 
      rr = (0:obj.conf.seg.ngrid)/obj.conf.seg.ngrid;
      obj.thres.xxgb = round(obj.mfi.sz(2)*rr); 
      obj.thres.yygb = round(obj.mfi.sz(1)*rr); 
      yy = obj.thres.yyg(2:end-1)-(obj.mfi.sz(1)+1)/2;
      xx = obj.thres.xxg(2:end-1)-(obj.mfi.sz(2)+1)/2;
      [YY,XX]=ndgrid(yy,xx);
      obj.thres.d2center = sqrt(XX.^2+YY.^2);
      obj.thres.isouter = obj.thres.d2center>=min(max(abs(yy)),max(abs(xx)));
      obj.thres.have_v = zeros(obj.conf.seg.ngrid);
      for ix=1:obj.conf.seg.ngrid
        for iy=1:obj.conf.seg.ngrid
          if mean(mean(obj.mfi.sm(obj.thres.yygb(iy)+1:obj.thres.yygb(iy+1),obj.thres.xxgb(ix)+1:obj.thres.xxgb(ix+1))))>0.1
            obj.thres.have_v(iy,ix) = 1;
          end
        end
      end
      r = obj.conf.seg.rneigh;
      q = uint8(ones(round(obj.mfi.sz/obj.conf.seg.ngrid)+2*r));
      q(1:round(size(q,1)/2),:)=3;
      obj.thres.norm_segcrit = max(vmlCalcSegCrit(q,3,[r+1 size(q,2)-r],[r+1 size(q,1)-r],r,0));
      c = [mean(obj.mfi.yy) mean(obj.mfi.xx)];
      d2sun = (obj.mfi.YY-c(1)).^2+(obj.mfi.XX-c(2)).^2;
      q = d2sun>obj.conf.r0sun_px^2 & d2sun<=obj.conf.rsun_px^2;
      q = double(q(any(q,2),any(q,1)));
      q = q+2*q.*(repmat(1:size(q,2),size(q,1),1)>size(q,2)/2);
      obj.thres.norm_segcrit_circ = max(vmlCalcSegCrit(uint8(q),3,[1 size(q,2)],[1 size(q,1)],r,0));
      
      %compute data for the SDF
      [~,xx,yy,obj.mfi.sdf.ksize]=vmlDownscale(double(obj.mfi.sm),obj.conf.sdf.sz);
      obj.mfi.sdf.yy = interp1(1:obj.mfi.sz(1),obj.mfi.yy,yy);
      obj.mfi.sdf.xx = interp1(1:obj.mfi.sz(2),obj.mfi.xx,xx);
      
      %clear sky irradiation
      tt = (obj.data.tmin:6/86400:obj.data.tmax)';
      [GHI, ~, DHI] = pvl_clearsky_ineichen(pvl_maketimestruct(...
        tt,obj.calib.UTC),obj.calib.model3D.Location);
      effDNI =  GHI-DHI;
      obj.data.IrrClear = [tt GHI effDNI DHI];
      
      %clear sky power model
      switch lower(plantid)
        case 'cavriglia'
          obj.data.PrelClear = [tt vml_prel_clear_cavriglia(obj,tt)];
      end
      
    end
    
    function print(obj,printlevel,txt)
      if obj.uidata.printlevel>=printlevel, disp(txt); end
    end
    
    function reset(obj)
      obj.calc = [];
      obj.calc.mvec = nan(2,length(obj.data.ti));
      obj.calc.thres_v = nan(obj.conf.seg.ngrid^2,length(obj.data.ti));
    end
    
    function x = imread(obj,j)
      %read a single frame
      while 1
        success = 1;
        try
          x = imread([obj.data.folder obj.data.files{j} obj.conf.img_ext]);
        catch
          success = 0;
        end
        if success, break; end
        obj.print(0,'error reading image, retrying after 5 seconds');
        pause(5);
      end
      if ~isempty(obj.oi.sky_area)
        x = x(obj.oi.sky_area(1):obj.oi.sky_area(2),obj.oi.sky_area(3):obj.oi.sky_area(4),:);
      end
    end
        
    function P1 = getP(obj,t)
      %get the power data for time(s) t
      if nargin<2, t = obj.data.ti;
      elseif any(t<obj.data.ti(1)-1), t=obj.data.ti(t); 
      end
      P1 = interp1(obj.data.P(:,1),obj.data.P(:,2),t);
    end
    
    function y = getIrr(obj,t)
      %get the irradiation data for time(s) t
      if nargin<2, t = obj.data.ti;
      elseif any(t<obj.data.ti(1)-1), t=obj.data.ti(t); 
      end
      y = interp1(obj.data.Irr(:,1),obj.data.Irr(:,2),t);
    end
    
    function y = getIrrClear(obj,t,col)
      %get the clear sky irradiation for time(s) t
      if nargin<2, t = obj.data.ti;
      elseif any(t<obj.data.ti(1)-1), t=obj.data.ti(t); 
      end
      if nargin<3, col=1; end
      y = interp1(obj.data.IrrClear(:,1),obj.data.IrrClear(:,1+col),t);
    end
    
    function y = getGHIClear(obj,t)
      if nargin<2, t = obj.data.ti; end
      y = obj.getIrrClear(t,1);
    end
    
    function y = getEffDNIClear(obj,t)
      if nargin<2, t = obj.data.ti; end
      y = obj.getIrrClear(t,2);
    end
    
    function y = getDHIClear(obj,t)
      if nargin<2, t = obj.data.ti; end
      y = obj.getIrrClear(t,3);
    end
    
    function y = getPrelClear(obj,t)
      %get the relative clear sky irradiation for time(s) t
      if nargin<2, t = obj.data.ti;
      elseif any(t<obj.data.ti(1)-1), t=obj.data.ti(t); 
      end
      y = interp1(obj.data.PrelClear(:,1),obj.data.PrelClear(:,2),t);
    end
    
    function y = getPersistErrP(obj,t,tpred)
      y = zeros(length(tpred),length(t));
      for i=1:length(tpred)
        t1 = t+tpred(i)/86400;
        y(i,:) = obj.getP(t).*obj.getPrelClear(t,tpred(i))./...
          obj.getPrelClear(t)-obj.getP(t1);
      end
    end
    
    function y = getPersistErrIrr(obj,t,tpred)
      if any(t<obj.data.ti(1)-1), t=obj.data.ti(t); end
      y = zeros(length(tpred),length(t));
      for i=1:length(tpred)
        t1 = t+tpred(i)/86400;
        y(i,:) = obj.getIrr(t).*obj.getIrrClear(t1)./...
          obj.getIrrClear(t)-obj.getIrr(t1);
      end
    end
    
    function j = getidx(obj,t,ttref)
      %get frame number(s) for given time(s)
      if nargin<3, ttref=obj.data.ti; end
      if ischar(t), t = floor(ttref(1))+time_str2num(t); 
      elseif t<ttref(1)-1, t=obj.data.ti(t);
      end
      [~,j] = min(abs(ttref-t));
    end
    
    function jj = find_daylight_range(obj,extend,thres)
      if nargin<2, extend=0; end
      if nargin<3, thres = obj.conf.daylight_thres; end
      jj = obj.data.IrrClear(:,2)>=max(obj.data.IrrClear(:,2))*thres;
      jstart = obj.getidx(obj.data.IrrClear(find(jj,1),1)-extend*max(obj.conf.pred.tpred)/86400);
      jend = obj.getidx(obj.data.IrrClear(find(jj,1,'last'),1));
      jj = jstart:jend;
    end
    
% end basic functions

%% cloud segmentation

    function varargout = seggrid2mfi(obj,varargin)
      [YYg,XXg] = ndgrid(obj.thres.yyg,obj.thres.xxg);
      [YY,XX] = ndgrid(1:length(obj.mfi.yy),1:length(obj.mfi.xx));
      VVg = XXg;
      varargout = cell(size(varargin));
      for i=1:numel(varargin)
        VVg(2:end-1,2:end-1) = varargin{i};
        VVg(1,:) = VVg(2,:); VVg(end,:) = VVg(end-1,:);
        VVg(:,1) = VVg(:,2); VVg(:,end) = VVg(:,end-1);
        y= interp2(XXg,YYg,VVg,XX,YY);
        y(~obj.mfi.sm) = NaN;
        varargout{i} = y;
      end
    end

    function loadframe(obj,j,redo_thres_calc)
      %load frame number j:
      % - read the image
      % - compute the median downscaled image
      % - identify opposite class neighbors
      % - optimize the cloud segmentation threshold map
      if nargin<3, redo_thres_calc=0; end
      if j<1 || j>length(obj.data.ti), obj.curseg=[]; return; end
      if isfield(obj.curseg,'j') && ~isempty(obj.curseg.j) && ...
          obj.curseg.j==j && ~redo_thres_calc
        return; 
      end
      obj.curseg.j = [];
      obj.curseg.x0 = obj.imread(j);
      obj.curseg.x1 = double(obj.curseg.x0);
      obj.curseg.thres.pre_cc = zeros(obj.mfi.sz);
      
      %detect clear sun by saturation and star shape
      sunp = obj.sunpos_im(j);
      [obj.curseg.clear_sun_flag,obj.calc.seg.yx_sun_detected(:,j),bestalpha_hex] = ...
        vmlDetectSun(obj.curseg.x0,obj.conf.sundetect,sunp);
      obj.calc.seg.clear_sun_flag(j) = obj.curseg.clear_sun_flag;

      if obj.curseg.clear_sun_flag == -1
        warning('sun could not be detected (outside image) ');
        obj.curseg.sm = obj.skymask_wo_sun(j,obj.mfi,obj.conf.r0sun_px);
      elseif obj.curseg.clear_sun_flag == -9
        warning('sun was incorrectly detected');
        obj.curseg.sm = obj.skymask_wo_sun(j,obj.mfi,obj.conf.rsun_px);
      elseif obj.curseg.clear_sun_flag == 1
        obj.curseg.sm = obj.skymask_wo_sun(j,obj.mfi,obj.conf.r0sun_px);
      else
        %correct bias close to sun
        sunp = obj.calc.seg.yx_sun_detected(:,j);
        d2sun = sqrt((obj.oi.YY-sunp(1)).^2+(obj.oi.XX-sunp(2)).^2);
        aa = (1:360)/180*pi;
        nflt = obj.conf.corr_sun_bias.nfilter;
        rr = obj.conf.r0sun_px-nflt+1:obj.conf.rsun_px+nflt;
        prgb=vmlImage2Polar(obj.oi,sunp,aa,rr,cat(3,obj.curseg.x1,obj.oi.sm));
        prgb = prgb(all(prgb(:,:,4)==1,2),:,1:3); %cut out directions which are not fully in sm
        pr2b = vmlRed2Blue(prgb);
        [~,jj] = sort(mean(pr2b(:,1+nflt:end-nflt),2));
        c0 = [0 0]; c_im = zeros([obj.oi.sz 2]); 
        c_im1 = zeros([size(pr2b) 2]);
        n = round(obj.conf.corr_sun_bias.clear_frac*size(pr2b,1));
        fcorr = obj.oi.sm; fcorr1 = ones(size(pr2b));
        for i=1:2
          c = filter(ones(1,1+2*nflt)/(1+2*nflt),1,mean(prgb(jj(1:n),:,i)));
          c = c(1+2*nflt:end);
          [c0(i),j1] = min(c); [cmax, j2] = max(c);
          if j2>j1, fcorr = zeros(size(fcorr)); break; end %no characteristic bias present
          c(j1+1:end) = c0(i); c(1:j2-1) = cmax;
          c_im(:,:,i) = interp1(0:obj.conf.rsun_px,...
            [cmax+zeros(1,1+obj.conf.r0sun_px) c],d2sun,'linear',c0(i));
          c_im1(:,:,i) = interp1(0:obj.conf.rsun_px,...
            [cmax+zeros(1,1+obj.conf.r0sun_px) c],repmat(rr,size(pr2b,1),1),'linear',c0(i));
          fcorr = min(fcorr,max(0,min(1,(obj.curseg.x1(:,:,i)-c0(i))./(c_im(:,:,i)-c0(i)))));
          fcorr1 = min(fcorr1,max(0,min(1,(prgb(:,:,i)-c0(i))./(c_im1(:,:,i)-c0(i)))));
        end
        for i=1:2
          obj.curseg.x1(:,:,i) = max(0,obj.curseg.x1(:,:,i))-fcorr.*(c_im(:,:,i)-c0(i));
        end
        if obj.curseg.clear_sun_flag>=3
          obj.curseg.sm = obj.skymask_wo_sun(j,obj.mfi,obj.conf.r0sun_px,sunp,bestalpha_hex);
        else
          obj.curseg.sm = obj.skymask_wo_sun(j,obj.mfi,obj.conf.r0sun_px,sunp);
        end
      end
      obj.curseg.x1 = uint8(obj.curseg.x1);
      
      sm = obj.curseg.sm;
      obj.curseg.x = vmlMedianDownscale(obj.curseg.x1,obj.conf.mfi.sz);
      obj.curseg.r2b = vmlRed2Blue(obj.curseg.x);
      obj.curseg.r2b(~sm) = NaN;
      
      %remove sunglare
      g2b = vmlGreen2Blue(obj.curseg.x);
      jj = find(sm & obj.curseg.x(:,:,2)<min(obj.curseg.x(:,:,1),obj.curseg.x(:,:,3)) & ...
        (obj.curseg.r2b>=0 | g2b<2*obj.curseg.r2b));
      n = min(length(jj),round(sum(sm(:))*obj.conf.seg.max_sunglare_remove));
      if n<length(jj)
        [~,jj1]=sort(obj.curseg.r2b(jj)-g2b(jj),1,'descend');
        jj = jj(jj1(1:n));
      end
      obj.curseg.r2b(jj) = NaN;
      
      %r2b map for current frame is now finalized
      r2b = obj.curseg.r2b;
      obj.curseg.d2sun = sqrt((obj.mfi.YY-sunp(1)).^2+(obj.mfi.XX-sunp(2)).^2);
      
      res = obj.conf.seg.r2b_resolution;
      resh = obj.conf.seg.r2b_hist_resolution;
      
      %identify clear sky areas close to sun, to be treated extra
      if obj.curseg.clear_sun_flag>=3
        jj = obj.curseg.d2sun<=obj.conf.rsun_px;
        r2b1 = obj.curseg.r2b; r2b1(~jj)=nan;
        r2b1 = r2b1(any(~isnan(r2b1),2),any(~isnan(r2b1),1));
        r2bmin = quantil(r2b1(:),obj.conf.seg.minmax_quantil);
        r2bmax = quantil(r2b1(:),1-obj.conf.seg.minmax_quantil);
        zzh = r2bmin-resh*1.5:resh:r2bmax+resh*2;
        r2b1q = median_filter_2d(vmlQuantize(r2b1,zzh),obj.conf.seg.rneigh,1,1);
        segcrit = vmlCalcSegCrit(r2b1q,length(zzh),[1 size(r2b1,2)],...
          [1 size(r2b1,1)],obj.conf.seg.rneigh,obj.conf.seg.halfgap)/obj.thres.norm_segcrit_circ;
        [s,j1] = max(segcrit);
        s = (s-obj.conf.seg.segcrit_close2sun(1))/diff(obj.conf.seg.segcrit_close2sun);
        if s>0 || obj.curseg.clear_sun_flag>=4
          if s>0
            obj.curseg.thres.pre_cc(jj) = -min(1,s)*(r2b(jj) < zzh(j1)); 
          end
          if obj.curseg.clear_sun_flag>=4
            if s>0, r=obj.conf.seg.rsun_fillclear_precc;
            else r=obj.conf.rsun_px; end
            obj.curseg.thres.pre_cc(obj.curseg.d2sun<=r)=-1;
          end
          jj = obj.curseg.d2sun<=obj.conf.seg.rsun_fillclear;
          r2b1 = r2b(obj.curseg.thres.pre_cc<0);
          r2blo = quantil(r2b1,obj.conf.seg.minmax_quantil);
          r2bhi = quantil(r2b1,1-obj.conf.seg.minmax_quantil);
          r2b(jj) = r2blo+(0:nnz(jj)-1)/(nnz(jj)-1)*(r2bhi-r2blo);
        end
      end
        
      ngrid = obj.conf.seg.ngrid;
      if all(isnan(obj.calc.thres_v(:,j))), redo_thres_calc = true; end
      if redo_thres_calc, obj.print(1,['cloud segmentation for frame #' num2str(j)]); end
      
      obj.curseg.thres.j0 = nan(ngrid);
      obj.curseg.thres.maxsegcrit = nan(ngrid);
      obj.curseg.thres.p_tanh_cumsum = nan(ngrid,ngrid,3);
      obj.curseg.thres.single_peak_crit = nan(ngrid);
      obj.curseg.thres.plot_data = cell(ngrid);
      obj.curseg.thres.isouter = obj.thres.isouter;
      obj.curseg.thres.v0_mask = nan(size(sm));
      obj.curseg.thres.v0lo = nan(ngrid);
      obj.curseg.thres.v0hi = nan(ngrid);
      obj.curseg.thres.hist_left_open_crit = nan(ngrid);
      obj.curseg.thres.hist_right_open_crit = nan(ngrid);
      
      obj.curseg.thres.norm_segcrit = obj.thres.norm_segcrit;
      obj.curseg.thres.yygb = obj.thres.yygb;
      obj.curseg.thres.xxgb = obj.thres.xxgb;
      
      %first identify r2b range to consider
      r2bmin = inf; r2bmax = -inf;
      for ix=1:ngrid
        for iy=1:ngrid
          r2b1 = r2b(obj.thres.yygb(iy)+1:obj.thres.yygb(iy+1),...
            obj.thres.xxgb(ix)+1:obj.thres.xxgb(ix+1));
          if mean(isnan(r2b1(:)))<=obj.conf.seg.tile_max_nan
            r2b1 = r2b1(~isnan(r2b1));
            obj.curseg.thres.v0lo(iy,ix) = quantil(r2b1,obj.conf.seg.minmax_quantil);
            obj.curseg.thres.v0hi(iy,ix) = quantil(r2b1,1-obj.conf.seg.minmax_quantil);
            r2bmin = min(r2bmin,obj.curseg.thres.v0lo(iy,ix)-obj.conf.seg.r2b_pad);
            r2bmax = max(r2bmax,obj.curseg.thres.v0hi(iy,ix)+obj.conf.seg.r2b_pad);
          end
        end
      end
      obj.curseg.thres.r2b_minmax = [r2bmin r2bmax];
      
      %initialize variables
      N = round((r2bmax-r2bmin)/res)+2;
      obj.curseg.thres.zz = r2bmin+(r2bmax-r2bmin-N*res)/2+res*(0:N);
      obj.curseg.thres.zzh = r2bmin-resh*1.5:resh:r2bmax+resh*2;
      obj.curseg.thres.ww = nan(N+1,ngrid,ngrid);
      obj.curseg.thres.X = tanh(bsxfun(@minus,obj.curseg.thres.zzh',obj.curseg.thres.zz)/res);
      %obj.curseg.thres.w_prior = vmlMakeSegPrior(obj.curseg.thres,obj.conf.seg);
      
      r2bq = median_filter_2d(vmlQuantize(r2b,obj.curseg.thres.zzh),obj.conf.seg.rneigh,1,1);
      obj.curseg.thres.r2bq = r2bq;
      obj.curseg.r2b_medianfilt = nan(size(r2bq));
      obj.curseg.r2b_medianfilt(r2bq>0) = obj.curseg.thres.zzh(r2bq(r2bq>0));
      
      ch1 = obj.curseg.x(:,:,2);
      ch1_filt = median_filter_2d(ch1,obj.conf.seg.rneigh,0,1);
      obj.curseg.thres.localvar = max(1,moving_avg_2d((ch1-ch1_filt).^2,obj.conf.seg.rneigh,1));
      
      if any(abs(obj.curseg.clear_sun_flag)==[1 2]), jj=true(obj.mfi.sz);
      else jj = obj.curseg.d2sun>obj.conf.rsun_px; end
      jj = find(jj(:));
      
      %pre_cc based on local variation
      llv = log(obj.curseg.thres.localvar(jj));
      th = obj.conf.seg.locvarthres;
      jj1 = llv<th(1);
      obj.curseg.thres.pre_cc(jj(jj1)) = -obj.conf.seg.locvarimpact(1)*(th(1)-llv(jj1))/th(1);
      jj1 = llv>th(1);
      obj.curseg.thres.pre_cc(jj(jj1)) = obj.conf.seg.locvarimpact(2)*min(1,(llv(jj1)-th(1))/diff(th));
      
      %pre_cc based on r2b prior
      r2b1 = r2b(jj);
      th = obj.conf.seg.prior_thres;
      jj1 = r2b1<th(2);
      obj.curseg.thres.pre_cc(jj(jj1)) = min(obj.curseg.thres.pre_cc(jj(jj1)),...
        -obj.conf.seg.prior_impact(1)*min(1,(th(2)-r2b1(jj1))/(th(2)-th(1))));
      jj1 = r2b1>th(3);
      obj.curseg.thres.pre_cc(jj(jj1)) = max(obj.curseg.thres.pre_cc(jj(jj1)),...
        obj.conf.seg.prior_impact(2)*min(1,(r2b1(jj1)-th(3))/(th(4)-th(3))));
      
      if redo_thres_calc
        for ix=1:ngrid
          for iy=1:ngrid
            jjy = obj.thres.yygb(iy)+1:obj.thres.yygb(iy+1);
            jjx = obj.thres.xxgb(ix)+1:obj.thres.xxgb(ix+1);
            d2sun1 = sqrt((obj.mfi.yy(obj.thres.yyg(iy+1))-sunp(1))^2+(obj.mfi.xx(obj.thres.xxg(ix+1))-sunp(2))^2);
            if d2sun1<=obj.conf.seg.rsun_maxouter, obj.curseg.thres.isouter(iy,ix) = 0; end
            
            if mean(mean(isnan(r2b(jjy,jjx))))>obj.conf.seg.tile_max_nan
              obj.curseg.thres.v0_mask(jjy,jjx) = -2;
            else
              obj.curseg.thres = vmlMakeSegHist(obj.curseg.thres,obj.conf.seg,ix,iy);
              
              if isnan(obj.curseg.thres.j0(iy,ix))
                obj.curseg.thres.v0_mask(jjy,jjx) = -1;
              else
                obj.curseg.thres.v0_mask(jjy,jjx) = obj.curseg.r2b(jjy,jjx)...
                  >=obj.curseg.thres.zzh(obj.curseg.thres.j0(iy,ix));
              end
            end
          end
        end
        
        v = vmlSmoothThresMap(obj.curseg.thres,obj.thres.have_v,obj.conf.seg);          
        obj.curseg.thres.v = extrapolate2nan(v);
        obj.calc.thres_v(:,j) = obj.curseg.thres.v(:);
      else
        obj.curseg.thres.v = reshape(obj.calc.thres_v(:,j),ngrid,ngrid);
      end
      
      [obj.curseg.thres.vmfi,obj.curseg.thres.vsingle_peak_crit] = ...
        obj.seggrid2mfi(obj.curseg.thres.v,obj.curseg.thres.single_peak_crit);

      obj.curseg.cc = double(obj.curseg.r2b>obj.curseg.thres.vmfi);
      obj.curseg.cc(isnan(obj.curseg.r2b) | ~sm) = NaN;
      if obj.curseg.clear_sun_flag>=4
        obj.curseg.cc(obj.curseg.d2sun<=obj.conf.seg.rsun_fillclear) = 0;
      end
      obj.curseg.sdf = vmlUpscale(vmlMakeSDF(obj.curseg.cc,obj.conf.sdf),obj.mfi.sz);
      obj.curseg.sdf(~obj.mfi.sm) = nan;
      
      obj.curseg.j = j;      
    end

% end cloud segmentation

%% update clear / cloudy to power / irr scaling

    function y = getCliP(obj,t)
      if nargin<2, t = obj.data.ti; end
      y = obj.getP(t)/(obj.getPrelClear(t)*obj.curscale.p_clear);
      y = 2*tanh(6*(y-obj.curscale.irr_th(1))/diff(obj.curscale.irr_th)-3)-1;
    end
    
    function y = getCliIrr(obj,t)
      if nargin<2, t = obj.data.ti; end
      y = obj.getIrr(t)/obj.getIrrClear(t);
      y = 2*tanh(6*(y-obj.curscale.irr_th(1))/diff(obj.curscale.irr_th)-3)-1;
    end

    function update_scale(obj)
      assert(obj.curseg.j==obj.curmot.j,'curseg and curmot need to be in sync');
      if isempty(obj.curscale)
        error('scaling was not initialized (probably startframe was not called)');
      end
      
      j = obj.curseg.j;
      irrel = obj.getIrr(j)/obj.getIrrClear(j);
      
      obj.calc.dhi_est(j) = nan;
      
      if irrel>=obj.conf.scale.irrelmin4upd
        if obj.curseg.clear_sun_flag==4 %clear
          obj.calc.dhi_est(j) = obj.getIrr(j)-obj.getEffDNIClear(j);
          
          obj.curscale.buf.irr_clear = [obj.curscale.buf.irr_clear irrel];
          if length(obj.curscale.buf.irr_clear)>=obj.conf.scale.buflen_min
            obj.curscale.buf.irr_clear(1:length(obj.curscale.buf.irr_clear)-...
              obj.conf.scale.buflen_max) = [];
            obj.curscale.irr_clear = median(obj.curscale.buf.irr_clear);
            if obj.curscale.irr_clear<obj.curscale.irr_cloudy
              warning('bad scaling: irr_clear < irr_cloudy');
              obj.curscale.irr_clear=obj.curscale.irr_cloudy;
            end
          end
          
          plant_cc = obj.calc_avgcc_plant(j,obj.curmot.sdf>0);
          if all(plant_cc<=obj.conf.scale.plantccmax4upd)
            %no cloud possibly on the plant: update power scale as well
            prel = obj.getP(j)/obj.getPrelClear(j)/irrel;
            obj.curscale.buf.p_clear = [obj.curscale.buf.p_clear prel];
            if length(obj.curscale.buf.p_clear)>=obj.conf.scale.buflen_min
              obj.curscale.buf.p_clear(1:length(obj.curscale.buf.p_clear)-...
                obj.conf.scale.buflen_max) = [];
              obj.curscale.p_clear = median(obj.curscale.buf.p_clear);
            end
          end
        elseif any(abs(obj.curseg.clear_sun_flag)==[1 2]) %cloudy
          obj.calc.dhi_est(j) = obj.getIrr(j);
          
          obj.curscale.buf.irr_cloudy = [obj.curscale.buf.irr_cloudy irrel];
          if length(obj.curscale.buf.irr_cloudy)>=obj.conf.scale.buflen_min
            obj.curscale.buf.irr_cloudy(1:length(obj.curscale.buf.irr_cloudy)-...
              obj.conf.scale.buflen_max) = [];
            obj.curscale.irr_cloudy = median(obj.curscale.buf.irr_cloudy);
            if obj.curscale.irr_clear<obj.curscale.irr_cloudy
              warning('bad scaling: irr_cloudy > irr_clear');
              obj.curscale.irr_cloudy=obj.curscale.irr_clear;
            end
          end
        end
      end
      for fn={'p_clear' 'irr_clear' 'irr_cloudy'}
        obj.calc.scale.(fn{1})(j) = obj.curscale.(fn{1});
      end
      
      obj.curscale.irr_th = [...
        min(obj.curscale.irr_clear,max(0.6,obj.curscale.irr_cloudy))-1e-8 ...
        max(obj.curscale.irr_cloudy,min(0.9,obj.curscale.irr_clear))+1e-8];
      obj.calc.cli.p(j) = obj.getCliP(j);
      obj.calc.cli.Irr(j) = obj.getCliIrr(j);
    end

% end update scaling

%% update predictions

    function calc_curpred(obj,j,q,i_t,clipred,cli2v,vpersist)
      w = exp(-obj.conf.pred.(q).sc2w(i_t,2)*obj.curpred.(q)(i_t).sc);
      w = w/sum(w);
      %record weights of the different expert groups
      [obj.calc.w_experts.(q){i_t}.cc2class(:,j),...
        obj.calc.w_experts.(q){i_t}.areas(:,j), ...
        obj.calc.w_experts.(q){i_t}.mvec_scale(:,j), ...
        obj.calc.w_experts.(q){i_t}.mvec(:,j), ...
        obj.calc.w_experts.(q){i_t}.shapechange(:,j)] = ...
        vmlSumsExceptDim(obj.curpred.p(i_t).w,obj.curpred.nexperts.(q));

      cli1 = (1+tanh(6*(w'*clipred)-3))/2;
      obj.curpred.(q)(i_t).cli = cli1;
      v1 = cli2v*(obj.curscale.irr_cloudy+cli*(obj.curscale.irr_clear-obj.curscale.irr_cloudy));
      obj.curpred.(q)(i_t).v = v1;
      vpred1 = [v1; vpersist];
      w = exp(-obj.conf.pred.(q).sc2w_agg(i_t,2)*obj.curpred.(q)(i_t).scagg);
      w = w/sum(w);
      obj.curpred.(q)(i_t).vagg = w'*vpred1;
      obj.calc.w_experts.(q){i_t}.agg(:,j) = w;
    end
    
    function update_pred(obj)
      j = obj.curmot.j;
      mv = obj.calc.mvec(:,max(1,j+1-obj.conf.pred.nn_mvec(end)):j);
      if any(isnan(mv(:))) || all(mv(:)==0), return; end
      pcur = obj.getP(j); irrcur = obj.getIrr(j); 
      irc0 = obj.getIrrClear(j); prc0 = obj.getPrelClear(j);
      
      for i_t=1:length(obj.conf.pred.tpred)
        tpred = obj.conf.pred.tpred(i_t);
        irc1 = obj.getIrrClear(obj.data.ti(j)+tpred/86400);
        prc1 = obj.getPrelClear(obj.data.ti(j)+tpred/86400);
        p_persist = pcur*prc1/prc0;
        irr_persist = irrcur*irc1/irc0;
        
        %compute experts' predictions
        ppred = obj.curpred.p(i_t).sc+nan;
        irrpred = obj.curpred.irr(i_t).sc+nan;
        k_p = 0; k_irr = 0; bb=[];
        for i_m = 1:length(obj.conf.pred.t_shapechange)
          for i_v = 1:length(obj.conf.pred.nn_mvec)
            for i_s = 1:length(obj.conf.pred.mvec_scale)
              n_mvec1 = obj.conf.pred.nn_mvec(i_v);
              mvec1 = obj.calc.mvec(:,j+1-n_mvec1:j);
              mvec1(:,all(mvec1==0,1)) = []; %remove zeros (typically uniform sky)
              sdf1 = obj.curmot.sdf+obj.curmot.m*min(tpred*...
                obj.conf.pred.f_shapechange(i_m),obj.conf.pred.t_shapechange(i_m));
              if ~isempty(mvec1)
                sdf1 = vmlShift(obj.mfi,mean(mvec1,2)*tpred*obj.conf.pred.mvec_scale(i_s),sdf1);
              end
              cc = sdf1>0;
              [p1,bb] = obj.calc_avgcc_plant(obj.data.ti(j)+tpred/86400,cc,[],bb);
              p1 = vml_cc2class(p1,obj.conf.pred.cc2class_cutoff);
              ppred(k_p+(1:length(p1))) = p1; k_p = k_p+length(p1);
              irr1 = vml_cc2class(obj.calc_avgcc_sun(obj.data.ti(j)+tpred/86400,cc),...
                obj.conf.pred.cc2class_cutoff);
              irrpred(k_irr+(1:length(irr1))) = irr1; k_irr = k_irr+length(irr1);
            end
          end
        end
        
        %aggregate experts' predictions to power / irr prediction
        obj.calc_curpred('p',j,i_t,ppred,prc1*obj.curscale.p_clear,p_persist);
        obj.calc_curpred('irr',j,i_t,irrpred,irc1,irr_persist);
        
        %store predictions in buffer for later weight update
        obj.curpred.p(i_t).buf = [obj.curpred.p(i_t).buf ppred];
        obj.curpred.irr(i_t).buf = [obj.curpred.irr(i_t).buf irrpred];
        obj.curpred.p(i_t).bufagg = [obj.curpred.p(i_t).bufagg ppred1];
        obj.curpred.irr(i_t).bufagg = [obj.curpred.irr(i_t).bufagg irrpred1];
        obj.curpred.tbuf{i_t} = [obj.curpred.tbuf{i_t} obj.data.ti(j)];
        
        %update scores
        while obj.data.ti(j)>=obj.curpred.tbuf{i_t}(1)+tpred/86400
          t1 = obj.curpred.tbuf{i_t}(1)+tpred/86400;
          if all(~isnan(obj.curpred.p(i_t).buf(:,1))) %power
            sc = obj.curscale.p_clear*obj.getPrelClear(t1);
            y = obj.getP(t1)/sc;
            x = obj.curpred.p(i_t).buf(:,1)/sc;
            grad = x*(obj.curpred.p(i_t).w'*x-y);
            obj.curpred.p(i_t).w = vmlSimplexProj(obj.curpred.p(i_t).w-obj.conf.pred.eta_p(i_t)*grad);
            x = obj.curpred.p(i_t).bufagg(:,1)/sc;
            grad = x*(obj.curpred.p(i_t).wagg'*x-y);
            obj.curpred.p(i_t).wagg = vmlSimplexProj(obj.curpred.p(i_t).wagg-obj.conf.pred.eta_p_agg(i_t)*grad);
          end
          obj.curpred.p(i_t).buf(:,1) = [];
          obj.curpred.p(i_t).bufagg(:,1) = [];
          
          if all(~isnan(obj.curpred.irr(i_t).buf(:,1))) %irr
            sc = obj.curscale.irr_clear*obj.getIrrClear(t1);
            y = obj.getIrr(t1)/sc;
            x = obj.curpred.irr(i_t).buf(:,1)/sc;
            grad = x*(obj.curpred.irr(i_t).w'*x-y);
            obj.curpred.irr(i_t).w = vmlSimplexProj(obj.curpred.irr(i_t).w-obj.conf.pred.eta_irr(i_t)*grad);
            x = obj.curpred.irr(i_t).bufagg(:,1)/sc;
            grad = x*(obj.curpred.irr(i_t).wagg'*x-y);
            obj.curpred.irr(i_t).wagg = vmlSimplexProj(obj.curpred.irr(i_t).wagg-obj.conf.pred.eta_irr_agg(i_t)*grad);
          end
          obj.curpred.irr(i_t).buf(:,1) = [];
          obj.curpred.irr(i_t).bufagg(:,1) = [];
          obj.curpred.tbuf{i_t}(1) = [];
        end

        %store predictions
        obj.calc.pred.p.v(i_t,j) = obj.curpred.p(i_t).v;
        obj.calc.pred.p.vagg(i_t,j) = obj.curpred.p(i_t).vagg;
        obj.calc.pred.irr.v(i_t,j) = obj.curpred.irr(i_t).v;
        obj.calc.pred.irr.vagg(i_t,j) = obj.curpred.irr(i_t).vagg;
      end
      
      if false
      %evaluate shifted trajectory errors
      n = length(obj.conf.pred.tshift); m = length(obj.conf.pred.traj);
      p_err = zeros(m,n); p_tshift = p_err;
      irr_err = zeros(m,n); irr_tshift = irr_err;
      tt1 = (0:obj.conf.pred.traj(end));
      ttraj = obj.data.ti(j)+tt1/86400;
      pred_ptraj = interp1([0 obj.conf.pred.tpred],[pcur obj.curpred.p_agg],0:obj.conf.pred.traj(end));
      pred_irrtraj = interp1([0 obj.conf.pred.tpred],[irrcur obj.curpred.irr_agg],0:obj.conf.pred.traj(end));
      persist_traj = obj.getIrrClear(ttraj)/obj.getIrrClear(j);
      ptraj = obj.getP(ttraj); irrtraj = obj.getIrr(ttraj);
      j1 = zeros(1,m);
      for i=1:m
        j1(i) = find(tt1==obj.conf.pred.traj(i),1);
        p_err(i,:) = mean(feval(obj.conf.pred.errfct,pred_ptraj(1:j1(i))-ptraj(1:j1(i))));
        irr_err(i,:) = mean(feval(obj.conf.pred.errfct,pred_irrtraj(1:j1(i))-irrtraj(1:j1(i))));
        obj.calc.trajpred.p.persist(i,j) = mean(feval(obj.conf.pred.errfct,pcur*persist_traj(1:j1(i))-ptraj(1:j1(i))));
        obj.calc.trajpred.irr.persist(i,j) = mean(feval(obj.conf.pred.errfct,irrcur*persist_traj(1:j1(i))-irrtraj(1:j1(i))));
        XX = [ones(j1(i),1) tt1(1:j1(i))'];
        b = XX\ptraj(1:j1(i))'; obj.calc.trajpred.p.trend(i,j) = b(2);
        b = XX\irrtraj(1:j1(i))'; obj.calc.trajpred.irr.trend(i,j) = b(2);
        b = XX\pred_ptraj(1:j1(i))'; obj.calc.trajpred.p.pred_trend(i,j) = b(2);
        b = XX\pred_irrtraj(1:j1(i))'; obj.calc.trajpred.irr.pred_trend(i,j) = b(2);
      end
      for tsh = 1:obj.conf.pred.tshift(end)
        itsh = find(tsh<=obj.conf.pred.tshift,1);
        for d=[-1 1]
          for i=1:m
            ptraj = obj.getP(ttraj-d*tsh/86400);
            irrtraj = obj.getIrr(ttraj-d*tsh/86400);
            p_err1 = mean(feval(obj.conf.pred.errfct,pred_ptraj(1:j1(i))-ptraj(1:j1(i))));
            irr_err1 = mean(feval(obj.conf.pred.errfct,pred_irrtraj(1:j1(i))-irrtraj(1:j1(i))));
            if p_err1<p_err(i,itsh), p_err(i,itsh:end)=p_err1; p_tshift(i,itsh:end)=d*tsh; end
            if irr_err1<irr_err(i,itsh), irr_err(i,itsh:end)=irr_err1; irr_tshift(i,itsh:end)=d*tsh; end
          end
        end
      end
      obj.calc.trajpred.p.err(:,j) = p_err(:);
      obj.calc.trajpred.p.tshift(:,j) = p_tshift(:);
      obj.calc.trajpred.irr.err(:,j) = irr_err(:);
      obj.calc.trajpred.irr.tshift(:,j) = irr_tshift(:);
      obj.print(1,['frame #' num2str(j) ': 5 min ahead error = ' ...
        num2str(obj.calc.trajpred.p.err(1,j)) ' vs. persistance = ' ...
        num2str(obj.calc.trajpred.p.persist(1,j))])
      end
    end

% end update predictions

%% process sequence

    function x = get_cur_xof(obj)
      x = vmlNormalizedRedChannel(obj.curseg) + 1j*obj.curseg.cc;
      x(~obj.curseg.sm) = NaN;
    end

    function startframe(obj, j, oldscale)
      obj.print(1,['starting at frame #' num2str(j)  ' / ' obj.data.day ' ' datestr(obj.data.ti(j),'HH:MM:SS')]);
      if nargin<3, oldscale=[]; end
      
      for fn={'p_clear' 'irr_clear' 'irr_cloudy'}
        if ~isempty(oldscale), obj.curscale.(fn{1})=oldscale.(fn{1});
        else obj.curscale.(fn{1})=obj.conf.scale.default.(fn{1}); end
        obj.curscale.buf.(fn{1}) = obj.curscale.(fn{1});
      end
      
      nexperts = [length(obj.conf.pred.cc2class_cutoff) ...
        length(obj.conf.plantproj.extensions)*length(obj.conf.plantproj.heights) ...
        length(obj.conf.pred.mvec_scale) ...
        length(obj.conf.pred.nn_mvec) ...
        length(obj.conf.pred.t_shapechange)];
      obj.curpred.nexperts.p = nexperts; nexperts_p = prod(nexperts);
      nexperts(2) = length(obj.conf.pred.rsun_px);
      obj.curpred.nexperts.irr = nexperts; nexperts_irr = prod(nexperts);
      for i=1:length(obj.conf.pred.tpred)
        obj.curpred.p(i).sc = zeros(nexperts_p,1)+.5;
        obj.curpred.irr(i).sc = zeros(nexperts_irr,1)+.5;
        obj.curpred.p(i).scagg = [.5;.5]; 
        obj.curpred.irr(i).scagg = [.5;.5];
        obj.curpred.p(i).buf = []; obj.curpred.p(i).bufagg = []; 
        obj.curpred.irr(i).buf = []; obj.curpred.irr(i).bufagg = [];
        obj.curpred.tbuf{i} = [];
      end
      
      obj.loadframe(j);
      obj.curmot.xofbuf = obj.get_cur_xof;
      obj.curmot.sdfbuf = obj.curseg.sdf;
      obj.curmot.j = j;
      obj.curmot.v = [0 0];
      obj.curmot.count = 1;
      obj.curmot.sdf = obj.curseg.sdf;
      obj.update_scale;
    end
    
    function nextframe(obj)
      assert(~isempty(obj.curmot.xofbuf),'nextframe called without startframe');
      j = obj.curmot.j+1;
      obj.print(1,' ');
      obj.print(1,['processing frame #' num2str(j) ' / ' obj.data.day ...
        ' ' datestr(obj.data.ti(j),'HH:MM:SS')]);
      obj.loadframe(j);
      obj.compute_next_mot;
      obj.update_scale;
      %obj.update_pred;
    end
    
    function process(obj,oldscale,jstart,jend)
      if nargin<2, oldscale = []; end
      if nargin<3, jstart = []; end
      if nargin<4, jend = []; end
      if isempty(jstart) || isempty(jend)
        jj = obj.find_daylight_range(1);
        if isempty(jstart), jstart = jj(1); end
        if isempty(jend), jend = jj(end); end
      end
      obj.curpred = [];
      obj.curpred.jend = jend;
      obj.startframe(jstart,oldscale);
      obj.resume(jend);
    end
    
    function resume(obj,jend)
      if nargin<2 
        if obj.curmot.j<obj.curpred.jend, jend = obj.curpred.jend;
        else jj = obj.find_daylight_range(1); jend = jj(end); end
      end
      obj.curpred.jend = jend;
      while obj.curmot.j<jend, obj.nextframe; end
    end
    
% end process sequence

%% plotting functionality
    
    function newfig(obj,figno,tit)
      if isempty(obj.uidata.varname)
        obj.uidata.varname = '<local>';
        for v = evalin('base','who')' 
          if strcmp(evalin('base',['class(' v{1} ')']),'vmlSeq');
            if (evalin('base',[v{1} '.uidata.id'])==obj.uidata.id)
              obj.uidata.varname = v{1};
              break;
            end
          end
        end
      end
      obj.uidata.open_figs = unique([obj.uidata.open_figs figno]);
      figure(figno); clf;
      set(gcf,'Name',[tit ' @ ' obj.uidata.varname]);%,'NumberTitle','off')
    end
    
    function plotSun(obj,j,col,marker,marker2)
      if nargin<3, col = 'k'; end
      if nargin<4, marker = []; end
      if nargin<5, marker2 = []; end
      tt = (0:.1:360)*pi/180;
      yx = sunpos_im(obj,j);
      xx = yx(2)+cos(tt)*obj.conf.rsun_px;
      xx(xx<1 | xx>obj.oi.sz(2)) = NaN;
      yy = yx(1)+sin(tt)*obj.conf.rsun_px;
      yy(yy<1 | yy>obj.oi.sz(1)) = NaN;
      if ~isempty(marker)
        line(yx(2),yx(1),'color',col,'marker',marker); 
      end
      if ~isempty(marker2) && isfield(obj.calc,'yx_sun_detected') && ...
          size(obj.calc.yx_sun_detected,2)>=j && all(~isnan(obj.calc.yx_sun_detected(:,j))) &&...
          all(obj.calc.yx_sun_detected(:,j)>0)
        line(obj.calc.yx_sun_detected(2,j),obj.calc.yx_sun_detected(1,j),'color',col,'marker',marker2); 
      end
      line(xx,yy,'color',col);
    end
    
    function showframe(obj,j,show_sun_or_skymask)
      %show the frame, part which does not belong to the sky is
      %greened out
      if nargin<2, j=obj.curseg.j; end
      if nargin<3, show_sun_or_skymask = 0; end
      if ~isscalar(show_sun_or_skymask), sm = show_sun_or_skymask;
      elseif show_sun_or_skymask>=2, sm = obj.skymask_wo_sun(j);
      else sm = obj.oi.sm; 
      end
      x = vmlColorify(obj.imread(j),~sm,2,64);
      image(x); axis off;
      if isscalar(show_sun_or_skymask) && show_sun_or_skymask==1
        obj.plotSun(j,'g','.','x');
      end
      title([datestr(obj.data.ti(j),'HH:MM:SS') ' (#' num2str(j) ')']);
    end
    
    function ShowEcalib(obj)
      %produces plot to visually check the external calibration
      if ~isfield(obj.calib,'kabsch') || isempty(obj.calib.kabsch)
        error('no external calibration was executed on this vmlSeq object');
      end
      figno = 1323; obj.newfig(figno,'External calibration');
      subplot(3,2,5:6);
      plot(obj.calib.kabsch.dt,obj.calib.kabsch.im_model(1,:),'-b.',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_model(2,:),'-r.',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_detect(1,:),'bo',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_detect(2,:),'ro');
      title(obj.data.day,'interpreter','none');
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
        figure(figno);
        subplot(3,2,1:4); 
        if ~show_processed
          obj.showframe(obj.calib.kabsch.frame_nos(j));
        else
          x = obj.im4sundetect(obj.calib.kabsch.frame_nos(j));
          imagesc(x); axis off;
          title([datestr(obj.data.ti(obj.calib.kabsch.frame_nos(j)),'HH:MM:SS') ' (#' num2str(obj.calib.kabsch.frame_nos(j)) ')']);
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
      plot(obj.data.P(:,1),obj.getP(obj.data.P(:,1)+toffset),'b',...
        obj.data.ti,obj.getP(obj.data.ti+toffset),'b.');
      ylim = get(gca,'ylim'); ylim(1) = 0;
      set(gca,'ylim',ylim);
      grid on;
      datetickzoom;
    end
    
    function plotirr(obj)
      %plot the power profiles
      plot(obj.data.ti,obj.getIrrClear(obj.data.ti),'r',...
        obj.data.ti,obj.getIrr(obj.data.ti),'b.');
      ylim = get(gca,'ylim'); ylim(1) = 0;
      set(gca,'ylim',ylim);
      grid on;
      datetickzoom;
    end
    
    function showtraj(obj,j)
      %plot the motion trajectory for frame j
      obj.showframe(j,1);
      mvec1 = obj.calc.mvec(:,j);
      if all(~isnan(mvec1))
        n = 300; n1 = 60;
        p = obj.sunpos_plane(j);
        pp = repmat(p,1,n)-mvec1*((1:n));
        yx = obj.camworld2im([pp;ones(1,n)]);
        hold on;
        plot(yx(2,:),yx(1,:),'r',yx(2,n1:n1:n),yx(1,n1:n1:n),'ro');
        hold off;
      end
      title([datestr(obj.data.ti(j),'HH:MM:SS') ' #' num2str(j)]);
    end
    
    function showcc(obj,j,cc)
      if nargin<3,
        j = obj.curmot.j;
        cc = obj.curmot.cc;
        cc(isinf(obj.curmot.age_cc)) = NaN;
      else
        obj.loadframe(j);
      end
      x = obj.curseg.x;
      x = vmlColorify(x,~obj.mfi.sm,1:3,-255);
      x = vmlColorify(x,isinf(obj.curmot.age_cc),1:3,40);
      image(obj.mfi.xx,obj.mfi.yy,x); axis off;
      title([datestr(obj.data.ti(j),'HH:MM:SS') ' (#' num2str(j) ')']);
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5]);
      hold on;
      contour(obj.mfi.XX,obj.mfi.YY,cc,[.5 .5],'g','linewidth',1);
      contour(obj.mfi.XX,obj.mfi.YY,isnan(cc),[.5 .5],'k','linewidth',1);
      hold off;
    end
    
    function plotplant(obj,j,h,ext)
      if nargin<3, h=1000; end
      if nargin<4, ext=0; end
      p = obj.plant2im(j,h,ext);
      line(p(2,:),p(1,:),'color','g');
    end
    
    function showplant(obj,j,h,ext)
      if nargin<3, h=1000; end
      if nargin<4, ext=0; end
      obj.showframe(j);
      obj.plotplant(j,h,ext);
    end
    
    function Show(obj, method)
      if nargin<2, method = 'showtraj'; end
      figno = 1321; obj.newfig(figno,'PV power and irradiation');
      ax(1) = subplot(2,4,3:4); obj.plotpow;
      title(obj.data.day,'interpreter','none');
      ax(2) = subplot(2,4,7:8); obj.plotirr;
      linkaxes(ax,'x');
      j = find(all(~isnan(obj.calc.mvec),1),1);
      if isempty(j), j=1; end
      show_upd(j);
      dcmobj = datacursormode(gcf); 
      set(dcmobj,'UpdateFcn',@show_upd);
      function out = show_upd( j, event_obj)
        if nargin>=2
          [~,j] = min(abs(event_obj.Position(1)-obj.data.ti));
        end
        figure(figno);
        subplot(2,4,[1 2 5 6]); 
        obj.(method)(j);
        drawnow;
        if nargin<2, out=[]; 
        else out = datestr(event_obj.Position(1),'HH:MM:SS'); end
      end
    end
    
    function ShowMvec(obj)
      if all(isnan(obj.calc.mvec(:))), error('no motion vectors available'); end
      figno = 1324; obj.newfig(figno,'Motion vectors');
      subplot(3,2,5:6);
      jj = find(all(~isnan(obj.calc.mvec)));
      if isempty(jj), error('no motion vectors computed yet'); end
      plot(obj.data.ti(jj),obj.calc.mvec(1,jj),'-b.',obj.data.ti(jj),obj.calc.mvec(2,jj),'-r.');
      datetickzoom;
      title(obj.data.day,'interpreter','none');
      show_mvec_upd(jj(1));
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@show_mvec_upd)
      function out = show_mvec_upd( j, event_obj)
        if nargin>=2, [~,j] = min(abs(event_obj.Position(1)-obj.data.ti)); end
        figure(figno);
        subplot(3,2,1:4); 
        obj.showtraj(j);
        drawnow;
        if nargin<2, out=[]; 
        else out = datestr(event_obj.Position(1),'HH:MM:SS'); end
      end
    end
    
    function showseg(obj,j)
      if nargin<2, 
        j = obj.curseg.j;
        assert(~isempty(j),'j not specified (previous computation was possibly aborted)');
      else
        obj.loadframe(j);
      end
      image(obj.mfi.xx,obj.mfi.yy,vmlColorify(obj.curseg.x,~obj.mfi.sm,1:3,0,1));
      axis off;
      axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5]);
      obj.plotSun(j);
      for p=(1:obj.conf.seg.ngrid-1)/obj.conf.seg.ngrid
        line(1+[p p]*(obj.oi.sz(2)-1),[1 obj.oi.sz(1)],'Color',[0 .5 0]);
        line([1 obj.oi.sz(2)],1+[p p]*(obj.oi.sz(1)-1),'Color',[0 .5 0]);
      end
      hold on;
      contour(obj.mfi.XX,obj.mfi.YY,obj.curseg.cc,[.5 .5],'g','linewidth',1);
      hold off;
      title([obj.data.day ', ' datestr(obj.data.ti(j),'HH:MM:SS') ' (#' num2str(j) ...
        ') / clear_sun_flag=' num2str(obj.curseg.clear_sun_flag)],'interpreter','none');
    end
    
    function ShowSeg(obj,j)
      if nargin<2, 
        j = obj.curseg.j;
        assert(~isempty(j),'j not specified (previous computation was possibly aborted)');
      else
        obj.loadframe(j);
      end
      figno = 1322; obj.newfig(figno,['Cloud thresholds, frame #' num2str(j) ' / ' datestr(obj.data.ti(j))]);
      cb_text = {'r2b','r2b mf','sm','glare','thres0','thres','th_sdf','thresv','spcrit','locvar','sdf','cont'};
      hcb = zeros(1,length(cb_text));
      htxt = uicontrol('Style','text','Position',[1 5 650 16]);
      for i = 1:length(hcb)
        hcb(i) = uicontrol('Style','checkbox','Value',i==length(hcb),...
          'Position',[(i-1)*47+1 24 47 16],'String',cb_text{i});
      end
      reset_zoom = 1;
      uicontrol('Style','pushbutton','Position',[5 45 60 20],'String','reset zoom','Callback',{@showSeg_upd1});
      auxplot_shown = [-1 -1];
      set(gcf,'toolbar','figure');
      set(hcb,'Callback',{@showSeg_upd}); %,hcb,hax
      showSeg_upd;

      function showSeg_click(varargin)
        p = get(gca,'currentpoint');
        [~,jx] = min(abs(p(1,1)-obj.mfi.xx));
        [~,jy] = min(abs(p(1,2)-obj.mfi.yy));
        ix = max(1,min(obj.conf.seg.ngrid,1+floor((p(1,1)-1)/(obj.oi.sz(2)-1)*obj.conf.seg.ngrid)));
        iy = max(1,min(obj.conf.seg.ngrid,1+floor((p(1,2)-1)/(obj.oi.sz(1)-1)*obj.conf.seg.ngrid)));
        if isfield(obj.curseg.thres,'single_peak_crit') && ~isnan((obj.curseg.thres.single_peak_crit(iy,ix)))
          strsingle_peak_crit = [', spcrit = ' num2str(obj.curseg.thres.single_peak_crit(iy,ix))];
        else
          strsingle_peak_crit = '';
        end
        set(htxt,'String',['xy = [' num2str(round(obj.mfi.xx(jx))) ', ' num2str(round(obj.mfi.yy(jy))) ...
           '] RGB = (' strjoin(cellfun(@num2str,squeeze(num2cell(obj.curseg.x(jy,jx,:))),'UniformOutput',0)',', ') ...
           ') R2B = ' num2str(obj.curseg.r2b(jy,jx),'%.2f') ' / ' ...
           'xytile = [' num2str(ix) ', ' num2str(iy) ']' strsingle_peak_crit ...
           ', localvar = ' num2str(obj.curseg.thres.localvar(jy,jx))]);
        if any([ix iy]~=auxplot_shown) && ~isempty(obj.curseg.thres.plot_data{iy,ix})
          zzh = obj.curseg.thres.zzh; 
          plot_data = obj.curseg.thres.plot_data{iy,ix};
          hh = plot_data{1}; segc = plot_data{2};
          wh = plot_data{3}; wc = plot_data{4}; w = plot_data{5};
          yyh = max(0,obj.curseg.thres.X*wh(2:end)+wh(1));
          yyc = max(0,obj.curseg.thres.X*wc(2:end)+wc(1));
          yy = obj.curseg.thres.X*w(2:end)+w(1);
          j0 = obj.curseg.thres.j0(iy,ix);
          figure(13281); clf; hold on;
          if obj.curseg.thres.single_peak_crit(iy,ix)>0, bar(zzh,hh,1,'m'); 
          else bar(zzh,hh,1,'r'); end
          bar(zzh,-segc,1,'b'); 
          plot(zzh,yyh,'r',zzh,-yyc,'b',zzh,yy,'k','linewidth',2);
          if ~isnan(j0), plot(zzh(j0),yy(j0),'ko','linewidth',2); end
          plot(obj.curseg.thres.v(iy,ix),interp1(zzh,yy,obj.curseg.thres.v(iy,ix)),'g*','linewidth',1);
          hold off; grid on;
        end
      end
      
      function showSeg_upd1(varargin)
        reset_zoom = 1;
        showSeg_upd;
      end
      
      function showSeg_upd(varargin)
        figure(figno);
        xlim = get(gca,'xlim'); ylim = get(gca,'ylim');
        if reset_zoom || (all(xlim==[0 1]) && all(ylim==[0 1]))
          xlim = [0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5];
          ylim = [0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5];
          reset_zoom = 0;
        end
        vv = get(hcb,'value');
        if vv{1} || vv{2}
          if vv{1}, x = obj.curseg.r2b; else x = obj.curseg.r2b_medianfilt; end
          x = max(0,x-quantil(x(:),0.05));
          x = repmat(uint8(min(1,x/quantil(x(:),0.95))*192),[1 1 3]);
        else
          x = obj.curseg.x;
        end
        x = vmlColorify(x,~obj.mfi.sm,1:3,128,1);
        if vv{3}
          x = vmlColorify(x,obj.mfi.sm & ~obj.curseg.sm,2,80);
        elseif vv{4}
          x = vmlColorify(x,isnan(obj.curseg.r2b) & obj.curseg.sm,2,80); 
        end
        if vv{5} && any(~isnan(obj.curseg.thres.v0_mask(:)))
          x = vmlColorify(x,obj.curseg.thres.v0_mask==-2,1:3,-80);
          x = vmlColorify(x,obj.curseg.thres.v0_mask==-1,1:3,40);
          x = vmlColorify(x,obj.curseg.thres.v0_mask==0,1,-40);
          x = vmlColorify(x,obj.curseg.thres.v0_mask==1,1,40);
        elseif vv{6}
          x = vmlColorify(x,obj.curseg.cc==0,1,-40);
          x = vmlColorify(x,obj.curseg.cc==1,1,40);
        elseif vv{7} && isfield(obj.curseg,'sdf')
          x = vmlColorify(x,obj.curseg.sdf<=0,1,-40);
          x = vmlColorify(x,obj.curseg.sdf>0,1,40);
        end
        ih = image(obj.mfi.xx,obj.mfi.yy,x); axis off;
        title([datestr(obj.data.ti(j),'HH:MM:SS') ' (#' num2str(j) ...
          ') / clear_sun_flag=' num2str(obj.curseg.clear_sun_flag)],'interpreter','none');
        axis([0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5 0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5]);
        set(ih,'ButtonDownFcn',@showSeg_click)
        obj.plotSun(j);
        for p=(1:obj.conf.seg.ngrid-1)/obj.conf.seg.ngrid
          line(1+[p p]*(obj.oi.sz(2)-1),[1 obj.oi.sz(1)],'Color',[0 .5 0]);
          line([1 obj.oi.sz(2)],1+[p p]*(obj.oi.sz(1)-1),'Color',[0 .5 0]);
        end
        for iy=1:obj.conf.seg.ngrid
          for ix=1:obj.conf.seg.ngrid
            if ~obj.thres.have_v(iy,ix)
              line(1+(ix-.5)/obj.conf.seg.ngrid*(obj.oi.sz(2)-1),1+(iy-.5)/obj.conf.seg.ngrid*(obj.oi.sz(1)-1),'Marker','x','Color',[0 .5 0])
            end
          end
        end
        ih = [];
        hold on;
        showthresval = (vv{8} && length(unique(obj.curseg.thres.vmfi(~isnan(obj.curseg.thres.vmfi))))>1);
        showsingle_peak_crit = vv{9};
        showlocalvar = vv{10};
        showsdf = vv{11};
        if showthresval || showsingle_peak_crit || showlocalvar || showsdf
          if showthresval, zz = obj.curseg.thres.vmfi; 
          elseif showsingle_peak_crit, zz = obj.curseg.thres.vsingle_peak_crit; 
          elseif showlocalvar, zz = log(0.5+obj.curseg.thres.localvar);
          elseif showsdf, zz = obj.curseg.sdf;
          end
          [~,ih] = contour(obj.mfi.XX,obj.mfi.YY,zz,20);
          set(gca,'clim',[min(zz(:)) max(zz(:))]);
          pos = get(gca,'Position'); 
          h = colorbar('vert');
          set(h,'AxisLocation','in');
          if ~showthresval && ~showsingle_peak_crit && showlocalvar
            set(h,'xticklabel',num2str(exp(get(h,'xtick')')-0.5,'%.1f'));
          end
          set(gca,'Position',pos);
        elseif vv{end} && length(unique(obj.curseg.cc(~isnan(obj.curseg.cc))))>1
          [~,ih] = contour(obj.mfi.XX,obj.mfi.YY,obj.curseg.cc,[.5 .5],'g','linewidth',1);
        end
        hold off;
        if ~isempty(ih), set(ih,'ButtonDownFcn',@showSeg_click); end
        set(gca,'xlim',xlim); set(gca,'ylim',ylim);
      end
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

    function p = sunpos_realworld(obj,t)
      %position of the sun in the world coordinates
      if t<obj.data.ti(1)-1, t = obj.data.ti(t); end
      t1 = [];
      [t1.year,t1.month,t1.day,t1.hour,t1.min,t1.sec]=datevec(t);
      t1.UTC = obj.calib.UTC;
      sun = sun_position(t1, obj.calib.model3D.Location);
      [x,y,z] = sph2cart(sun.azimuth*pi/180,(90-sun.zenith)*pi/180,1);
      p = [x y z]';
    end
    
    function y = cos_sun2panel(obj,t)
      [x,y,z] = sph2cart(obj.conf.panels.azimuth*pi/180,(90-obj.conf.panels.zenith)*pi/180,1);
      y = [x y z]*obj.sunpos_realworld(t);
    end
    
    function yx = sunpos_im(obj,t)
      %position of the sun in image pixels
      yx = obj.camworld2im(obj.calib.Rext'*obj.sunpos_realworld(t));
    end
    
    function p = sunpos_plane(obj,t)
      %position of the sun in the unit plane
      p = obj.calib.Rext'*obj.sunpos_realworld(t);
      p = p(1:2)/p(3);
    end
    
% end position of the sun
    
%% plant projection and averages

    function [p,outer_normal] = plant2plane(obj,t,h)
      % project plant shape onto sky image
      sunp = obj.sunpos_realworld(t); % line slope for projection
      %m(1:2)=obj.adjust_sun_pos_v2(m(1:2));
      if size(obj.conf.plant_coords,2)<1
        return;
      end
      p = obj.calib.Rext'*(bsxfun(@plus,...
        obj.conf.plant_coords,sunp*h/sunp(3)));
      p = bsxfun(@rdivide,p(1:2,:),p(3,:));
      p = p(:,[1:end 1]);
      dp = diff(p,[],2);
      sig = sign(sum(dp(1,:).*(p(2,1:end-1)+p(2,2:end))));
      dp = bsxfun(@rdivide,dp,sqrt(sum(dp.^2)));
      n0 = sig*[-dp(2,:);dp(1,:)];
      outer_normal = (n0+n0(:,[end 1:end-1]))/2;
      outer_normal = bsxfun(@rdivide,outer_normal,sqrt(sum(outer_normal.^2)));
      outer_normal = bsxfun(@times,outer_normal,1./dot(n0,outer_normal));
      outer_normal = outer_normal(:,[1:end 1]);
    end
    
    function p = plant2im(obj,t,h,ext)
      if nargin<4, ext=0; end
      [p,np] = obj.plant2plane(t,h);
      p = obj.camworld2im([p+np*ext;ones(1,size(p,2))]);
    end
    
    function y = calc_avgcc_sun(obj,t,cc)
      sunp = obj.sunpos_im(t);
      d2sun = sqrt((obj.mfi.YY-sunp(1)).^2+(obj.mfi.XX-sunp(2)).^2);
      y = zeros(length(obj.conf.pred.rsun_px),1);
      for i=1:length(y)
        y(i) = mean(cc(d2sun<=obj.conf.pred.rsun_px(i)));
      end
    end
    
    function [y, bb] = calc_avgcc_plant(obj,t,cc,ee,bb)
      if nargin<4, ee=0; end
      if nargin<5, bb={}; end
      if isempty(ee), ee = obj.conf.plantproj.extensions; end
      has_bb = ~isempty(bb);
      hh = obj.conf.plantproj.heights;
      y = zeros(length(ee)*length(hh),1);
      k = 0;
      j1 = ~isnan(obj.mfi.ppx(:,1));
      ccv = cc(j1);
      for ih=1:length(hh)
        if ~has_bb, [p0,n0p] = obj.plant2plane(t,hh(ih)); end
        for i=1:length(ee)
          k = k+1;
          if has_bb, y(k) = mean(ccv(bb{k})); else
            p = p0+n0p*ee(i);
            b=inpolygon(obj.mfi.ppx(j1,1),obj.mfi.ppx(j1,2),p(1,:),p(2,:));
            y(k) = mean(ccv(b));
            bb{k} = b;
          end
        end
      end
    end

% end plant projection and averages
    
%% sky masks
    
    function sm = skymask_wo_sun(obj,j,iminfo,rsun_px,yx_sunpos,bestalpha_hex)
      %sky mask without sun, frame j (taking out a circle around 
      %the sun position)
      if nargin<3, iminfo = obj.oi; end
      if nargin<4, rsun_px = obj.conf.rsun_px; end
      if nargin<5, yx_sunpos = sunpos_im(obj,j); end
      if nargin<6, bestalpha_hex = []; end
      sm = iminfo.sm;
      XX = iminfo.XX-yx_sunpos(2); YY = iminfo.YY-yx_sunpos(1);
      sm(XX.^2+YY.^2<=rsun_px^2) = false;
      dx = iminfo.xx(2)-iminfo.xx(1);
      if ~isempty(bestalpha_hex)
        for a=0:60:120
          sm(XX.^2+YY.^2<=obj.conf.rsun_px^2 & ...
            abs(XX*sind(bestalpha_hex+a)-YY*cosd(bestalpha_hex+a))<dx) = false;
        end
      end
    end
    
    function w = weights4oflow(obj,j,mvec1,iminfo)
      %compute sky mask for the part of the sky where the clouds
      %come from
      if nargin<4, iminfo = obj.mfi; end
      
      w = [];
      if all(mvec1==0) || obj.conf.oflow.low_weight>=1, return; end
      
      w = ones(iminfo.sz);
      mvec1 = mvec1(:)/sqrt(sum(mvec1.^2));
      yx = sunpos_plane(obj,j);
      mvec1rot = [-mvec1(2);mvec1(1)];
      llim = yx'*mvec1rot-obj.conf.oflow.heigh_weight_width; 
      hlim = yx'*mvec1rot+obj.conf.oflow.heigh_weight_width;
      if obj.conf.oflow.heigh_weight_cover_center
        if llim<0 && hlim<0, hlim=0; elseif llim>0 && hlim>0, llim=0; end
      end
      ppxproj = iminfo.ppx*mvec1rot;
      w(ppxproj<llim | ppxproj>hlim) = obj.conf.oflow.low_weight;
      w(iminfo.ppx*mvec1-yx'*mvec1>0) = obj.conf.oflow.low_weight;
    end
      
% end sky masks
    
%% sun detection / external calibration (Kabsch algorithm)
    
    function x = im4sundetect(obj,j)
      %compute image for automatic detection of the sun
      x = sum(obj.imread(j),3)/(3*255);
      x(~(obj.oi.sm | obj.oi.sm_border)) = 0.5;
      x = vmlSobelAbs(x);
      x(~obj.oi.sm) = 0;
      x=vmlGaussianBlur(x,obj.conf.sundetect.blur_px);
    end
    
    function yx = detect_saturated_sun(obj,j)
      %automatic detection of the (saturated) sun
      x = obj.im4sundetect(j);
      [~,jmax] = max(x(:));
      [p,q] = ind2sub(obj.oi.sz,jmax);
      yx = [p;q];
    end
    
    function R = kabsch(obj,R,frame_nos)
      %Kabsch algorithm for compution rotation matrix such that
      %detected sun position and acutal sun position coincide
      if nargin<2, R = []; end
      if nargin<3 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.data.ti); 
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
      obj.calib.kabsch.dt = obj.data.ti(frame_nos);
      obj.calib.kabsch.world_detect = PP;
      obj.calib.kabsch.realword = QQ;
      obj.calib.kabsch.im_detect = pp;
      obj.calib.kabsch.im_model = obj.camworld2im(R'*QQ);
      obj.calib.kabsch.R = R;
      obj.print(1,['#frames = ' num2str(length(frame_nos)) ', pixel RMSE = ' ...
        num2str(sqrt(mean(sum((pp-obj.calib.kabsch.im_model).^2,1))),'%.1f') ...
        ', RMSE on 1 m unit sphere = ' ...
        num2str(1000*sqrt(mean(sum((R*PP-QQ).^2,1))),'%.1f') ' mm']);
    end
    
% end sun detection

%% motion prediction and optical flow

    function compute_next_mot(obj)
      curmot1 = obj.curmot;
      j = curmot1.j+1;
      curmot1.j = j;
      curmot1.count = curmot1.count+1;
      if obj.curseg.j~=j, obj.loadframe(j); end

      %put last picture for optical flow and sdf into buffers
      curmot1.xofbuf = cat(3,curmot1.xofbuf,obj.get_cur_xof);
      if size(curmot1.xofbuf,3)>obj.conf.oflow.nframes
        curmot1.xofbuf = curmot1.xofbuf(:,:,end-obj.conf.oflow.nframes+1:end);
      end
      curmot1.sdfbuf = cat(3,curmot1.sdfbuf,obj.curseg.sdf);
      if size(curmot1.sdfbuf,3)>obj.conf.sdf.nframes
        curmot1.sdfbuf = curmot1.sdfbuf(:,:,end-obj.conf.sdf.nframes+1:end);
      end
            
      %compute optical flow
      if all(~isnan(obj.calc.mvec(:,j))) 
        curmot1.v = obj.calc.mvec(:,j);
      else
        obj.print(1,['optical flow for frame #' num2str(j)]);

        n = size(curmot1.xofbuf,3);
        tt = (obj.data.ti(j-n+1:j)-obj.data.ti(j))*86400;
        while any(abs(tt-(tt(1)+mean(diff(tt))*(0:length(tt)-1)))>obj.conf.oflow.t_noise_lim)
          n = n-1;
          tt = tt(2:end);
        end
        
        xof = curmot1.xofbuf(:,:,end+1-n:end);
        cc1 = imag(xof(~isnan(xof))); 
        if all(cc1==0) || all(cc1==1)
          v = [0 0]; %all clear or covered ==> no motion
        else
          v = curmot1.v;
          if curmot1.count<=2
            v = vmlOpticalFlow(obj.mfi,obj.conf.oflow,xof,tt,v);
          end
          w = obj.weights4oflow(j,v);
          v = vmlOpticalFlow(obj.mfi,obj.conf.oflow,xof,tt,v,w);
        end
        curmot1.v = v; obj.calc.mvec(:,j) = v(:);
      end
      
      %shift past data according to optical flow
      [curmot1.sdfbuf(:,:,1:end-1)] = vmlShift(obj.mfi,curmot1.v*...
        (obj.data.ti(j-1)-obj.data.ti(j))*86400,curmot1.sdfbuf(:,:,1:end-1));
      
      %cloud shape tracking
      if size(curmot1.sdfbuf,3)<obj.conf.sdf.nframes_min
        curmot1.sdf = obj.curseg.sdf;
        curmot1.m = obj.curseg.sdf*0;
      else
        obj.print(1,['shape tracking for frame #' num2str(j)]);
        jj = sum(~isnan(curmot1.sdfbuf),3)>=obj.conf.sdf.nframes_min;
        sz = size(curmot1.sdfbuf);
        sdf1 = reshape(curmot1.sdfbuf,sz(1)*sz(2),sz(3));
        tt = (obj.data.ti(j-sz(3)+1:j)-obj.data.ti(j))*86400;
        coeff=vmlPointwiseL1Regression(tt,sdf1(jj,:)',obj.conf.sdf.alpha,obj.conf.sdf.nframes_min);
        curmot1.sdf = nan(obj.mfi.sz); curmot1.sdf(jj)=coeff(1,:);
        curmot1.m = nan(obj.mfi.sz); curmot1.m(jj)=coeff(2,:);
      end
            
      obj.curmot = curmot1;
    end
  
% end motion prediction
  
  end 
end
