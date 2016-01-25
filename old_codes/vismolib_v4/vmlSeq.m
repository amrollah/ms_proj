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
    curpred;        %current prediction
  end
  
  methods
  
%% basic functions: constructor, imread, getP, ...

    function obj = vmlSeq(folder,hours,conf)
      if nargin<2, hours = []; end
      if nargin<3, conf = evalin('base','VMLCONF'); end
      
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
      obj.oi.sm_border = false(obj.oi.sz);
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
      plane_px1 = reshape(obj.oi.ppx,[obj.oi.sz 2]);
      d1 = plane_px1(2:end,:,1)-plane_px1(1:end-1,:,1);
      d2 = plane_px1(:,2:end,2)-plane_px1(:,1:end-1,2);
      obj.oi.dplane = min(min(abs(d1(~isnan(d1)))),min(abs(d2(~isnan(d2)))));
      assert(obj.oi.dplane >= 1/max(height,width));
      
      %load irradiation and power data
      if exist([conf.tmpfolder folder '\data.mat'],'file')
        load([conf.tmpfolder folder '\data.mat']);
      else
        run([conf.datafolder 'loaddata.m']);
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
      
      %assign title and day in date time format
      obj.uidata.strtitle = folder;
      obj.data.dt_day = datenum(obj.uidata.strtitle,'yyyy_mm_dd');
         
      %initialize data for motion vectors and power prediction
      obj.calc.mvec = nan(2,length(obj.data.ti));
      obj.calc.thres_v = nan(obj.conf.seg.ngrid^2,length(obj.data.ti));
      
      %do a median downscale of the image and compute the sky mask...
      [obj.mfi.sm,obj.mfi.xx,obj.mfi.yy] = ...
        vmlMedianDownscale(obj.oi.sm,obj.conf.mfi.nfilter);
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
      obj.thres.have_v = zeros(obj.conf.seg.ngrid);
      for ix=1:obj.conf.seg.ngrid
        for iy=1:obj.conf.seg.ngrid
          if any(any(obj.mfi.sm(obj.thres.yyg(iy):obj.thres.yyg(iy+2),obj.thres.xxg(ix):obj.thres.xxg(ix+2))))
            obj.thres.have_v(iy,ix) = 1;
          end
        end
      end
      rr = (0:obj.conf.seg.ngrid)/obj.conf.seg.ngrid;
      obj.thres.xxgb = round(obj.mfi.sz(2)*rr); 
      obj.thres.yygb = round(obj.mfi.sz(1)*rr); 
      yy = obj.thres.yyg(2:end-1)-(obj.mfi.sz(1)+1)/2;
      xx = obj.thres.xxg(2:end-1)-(obj.mfi.sz(2)+1)/2;
      [YY,XX]=ndgrid(yy,xx);
      obj.thres.isouter = XX.^2+YY.^2>=min(max(abs(yy)),max(abs(xx)))^2;
      
      %clear sky irradiation
      tt = (obj.data.tmin:6/86400:obj.data.tmax)';
      obj.data.IrrClear = [tt pvl_clearsky_ineichen(pvl_maketimestruct(tt, ...
        obj.calib.UTC),obj.calib.model3D.Location)];
    end
    
    function print(obj,txt,printlevel)
      if obj.uidata.printlevel>=printlevel, disp(txt); end
    end
    
    function x = imread(obj,j)
      %read a single frame
      x = imread([obj.data.folder obj.data.files{j} obj.conf.img_ext]);
      if ~isempty(obj.oi.sky_area)
        x = x(obj.oi.sky_area(1):obj.oi.sky_area(2),obj.oi.sky_area(3):obj.oi.sky_area(4),:);
      end
    end
        
    function P1 = getP(obj,t)
      %get the power data for time(s) t
      if any(t<obj.data.ti(1)-1), t=obj.data.ti(t); end
      P1 = interp1(obj.data.P(:,1),obj.data.P(:,2),t);
    end
    
    function y = getIrr(obj,t)
      %get the irradiation data for time(s) t
      if any(t<obj.data.ti(1)-1), t=obj.data.ti(t); end
      y = interp1(obj.data.Irr(:,1),obj.data.Irr(:,2),t);
    end
    
    function y = getIrrClear(obj,t)
      %get the clear sky irradiation for time(s) t
      if any(t<obj.data.ti(1)-1), t=obj.data.ti(t); end
      y = interp1(obj.data.IrrClear(:,1),obj.data.IrrClear(:,2),t);
    end
    
    function P1 = getPclearsky(obj,t)
      %compute the clear sky power model for time(s) t
      if any(t<obj.data.ti(1)-1), t=obj.data.ti(t); end
      tt = t(:)'-floor(obj.data.ti(1));
      P1 = sim_nn_nf(obj.Pclearsky.w,tt);
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
      for fn=fieldnames(obj.calc)'
        if ~isempty(obj.calc.(fn{1})) && any(~isnan(obj.calc.(fn{1})(:)))
          data.(fn{1}) = obj.calc.(fn{1});
        end
      end
    end
    
    function import(obj,data)
      for fn=fieldnames(data)', obj.calc.(fn{1}) = data.(fn{1}); end
    end
    
    function j = getidx(obj,t,ttref)
      %get frame number(s) for given time(s)
      if nargin<3, ttref=obj.data.ti; end
      if ischar(t), t = floor(ttref(1))+time_str2num(t); 
      elseif t<ttref(1)-1, t=obj.data.ti(t);
      end
      [~,j] = min(abs(ttref-t));
    end
    
% end basic functions

%% cloud segmentation
    
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
      obj.curseg.x = vmlMedianDownscale(obj.curseg.x0,obj.conf.mfi.nfilter);
      obj.curseg.r2b = vmlRed2Blue(obj.curseg.x);
      
      sunp = obj.sunpos_im(j);
      sunp_mfi = round(sunp.*(obj.mfi.sz./obj.oi.sz)');
      yx = round(sunp);
      if any(yx<1) || any(yx'>obj.oi.sz)
        obj.calc.satsun(j) = 0;
        obj.calc.star_crit(j) = 0; 
      else
        r = obj.conf.sundetect.r_roi_px+obj.conf.sundetect.r2_px;
        yy = max(1,yx(1)-r):min(obj.oi.sz(1),yx(1)+r);
        xx = max(1,yx(2)-r):min(obj.oi.sz(2),yx(2)+r);
        rsq = bsxfun(@plus,(yy-yx(1))'.^2,(xx-yx(2)).^2);
        x0 = mean(obj.curseg.x0(:,:,2:3),3);
        x1 = x0(yy,xx);
        x1 = vmlGaussianBlur(x1,obj.conf.sundetect.blur2_px);
        x2 = x1; x2(rsq>obj.conf.sundetect.r_roi_px^2) = inf;
        v_in = min(x2(:));
        [my,mx]=find(x1==v_in);
        my = round(mean(my)); mx = round(mean(mx));
        obj.calc.yx_sun_detected(:,j) = [yy(my);xx(mx)];
        obj.calc.yx_mfi_sun_detected(:,j) = round([yy(my);xx(mx)].*(obj.mfi.sz./obj.oi.sz)');
        obj.calc.sunpos_err(j) = norm(obj.calc.yx_sun_detected(:,j)-sunp);
        rsq=bsxfun(@plus,(yy'-yy(my)).^2,(xx-xx(mx)).^2);
        v_out = median(x1(rsq>=obj.conf.sundetect.r1_px^2 & rsq<=obj.conf.sundetect.r2_px^2));
        obj.calc.satsun(j) = v_out-v_in;
        star_rmax = min([obj.conf.sundetect.star_r2_px obj.oi.sz(1)-yy(my) yy(my)-1 obj.oi.sz(2)-xx(mx) xx(mx)-1]);
        if star_rmax<obj.conf.sundetect.star_r1_px || ...
          obj.calc.satsun(j)<obj.conf.sundetect.satsun_thres
          obj.calc.star_crit(j) = 0;
        else
          aa = (1:360)'/180*pi;
          rr = obj.conf.sundetect.star_r1_px:star_rmax;
          ww = (10.^-((rr-rr(1))'/(rr(end)-rr(1)))); ww = ww/sum(ww);
          zz=interp2(obj.oi.XX-xx(mx),obj.oi.YY-yy(my),mean(obj.curseg.x0(:,:,1:2),3),...
            bsxfun(@times,cos(aa),rr),bsxfun(@times,sin(aa),rr));
          zz = bsxfun(@minus,zz,mean(zz)); zz = bsxfun(@rdivide,zz,std(zz));
          obj.calc.star_crit(j) = 1-vmlHexFit(aa,ww,zz);
        end
      end
      
      if obj.calc.satsun(j)>=obj.conf.sundetect.satsun_thres
        rsun_px = obj.conf.rsun_px;
        if obj.calc.sunpos_err(j)>obj.conf.sundetect.maxerr_px
          obj.curseg.clear_sun_flag = -1;
        elseif obj.calc.star_crit(j)>=obj.conf.sundetect.star_thres_cov
          obj.curseg.clear_sun_flag = 3+(obj.calc.star_crit(j)>=obj.conf.sundetect.star_thres_clear);
        else
          obj.curseg.clear_sun_flag = 2;
          rsun_px = obj.conf.rsun_px_close;
        end
        obj.curseg.sm = obj.skymask_wo_sun(j,obj.mfi,rsun_px);
      else 
        obj.curseg.clear_sun_flag = 1;
        obj.curseg.sm = obj.mfi.sm;
      end
      obj.calc.clear_sun_flag(j) = obj.curseg.clear_sun_flag;
      sm = obj.curseg.sm;
      
      %remove sunglare
      g2b = vmlGreen2Blue(obj.curseg.x);
      center = obj.mfi.sz'/2;
%       sunp = obj.calc.yx_mfi_sun_detected(:,j);
      glare_ct = center + 0.7*(center-sunp_mfi);
      scale = 145/norm(center-sunp_mfi);
      [xx_g,yy_g] = ndgrid((1:obj.mfi.sz(1))-glare_ct(1),(1:obj.mfi.sz(2))-glare_ct(2));
      glare_mask = (xx_g.^2 + yy_g.^2) < (scale*1.2*obj.conf.sundetect.r_mfi_px)^2;
      
      glare_ct2 = sunp_mfi + 0.4*(sunp_mfi-center);
      [xx_g,yy_g] = ndgrid((1:obj.mfi.sz(1))-glare_ct2(1),(1:obj.mfi.sz(2))-glare_ct2(2));
      glare_mask2 = (xx_g.^2 + yy_g.^2) < (scale*.8*obj.conf.sundetect.r_mfi_px)^2;
      
%       glare_ct3 = sunp - 0.4*(sunp-center);
%       [xx_g,yy_g] = ndgrid((1:obj.mfi.sz(1))-glare_ct3(1),(1:obj.mfi.sz(2))-glare_ct3(2));
%       glare_mask3 = (xx_g.^2 + yy_g.^2) < (scale*.8*obj.conf.sundetect.r_mfi_px)^2;
      
      
      dir_vect =(center-sunp_mfi)./norm(center-sunp_mfi);
      m = [-dir_vect(2); dir_vect(1)];
      
      mid_st = sunp_mfi + dir_vect*obj.conf.sundetect.r_mfi_px;
      mid_end = glare_ct+ dir_vect*1.8*obj.conf.sundetect.r_mfi_px;
      p1 = mid_st - obj.conf.glare_rect_w * m;
      p2 = mid_st + obj.conf.glare_rect_w * m;      
      p3 = mid_end + obj.conf.glare_rect_w * m;      
      p4 = mid_end - obj.conf.glare_rect_w * m;
      
      P=round([p1,p2,p3,p4,p1]);
      [xx,yy] = ndgrid(1:obj.mfi.sz(1),1:obj.mfi.sz(2));
      glare_rect = inpolygon(xx,yy,P(1,:),P(2,:));
      
      [xx,yy] = ndgrid((1:obj.mfi.sz(1))-center(1),(1:obj.mfi.sz(2))-center(2));
      center_mask = (xx.^2 + yy.^2) < (6*obj.conf.sundetect.r_mfi_px)^2;

      obj.curseg.x_LUV = colorspace(['RGB->','Luv'], double(obj.curseg.x));
      luv_glare = obj.curseg.x_LUV(:,:,2)./obj.curseg.x_LUV(:,:,3);
      
%       t1=(obj.curseg.r2b<0 & luv_glare<0);
      t2=(obj.curseg.r2b>0 & (g2b<0));
      glare = (sm & (glare_mask | glare_mask2 | glare_rect) & (t2));
      jj=find(glare);
      
      obj.curseg.glare = glare;
      n = min(length(jj),round(sum(sm(:))*obj.conf.seg.max_sunglare_remove));
      if n<length(jj)
        [~,jj1]=sort(obj.curseg.r2b(jj),1,'descend');
        jj = jj(jj1(1:n));
      end
%       xflat = reshape(obj.curseg.x,size(obj.curseg.x,1)*size(obj.curseg.x,2),size(obj.curseg.x,3));
%       sky = (sm & center_mask & obj.curseg.r2b<0 & ~glare)>0;
%       obj.curseg.r2b(jj) = repmat(sky_r2b_avg,[numel(jj),1]);
      obj.curseg.r2b(jj) = nan;
%       xflat(jj,:) = repmat(sky_avg,[numel(jj),1]);
%       obj.curseg.x = reshape(xflat,size(obj.curseg.x,1),size(obj.curseg.x,2),size(obj.curseg.x,3));
            
      ngrid = obj.conf.seg.ngrid;
      if ~redo_thres_calc && any(~isnan(obj.calc.thres_v(:,j)))
        obj.curseg.thres = [];
        obj.curseg.thres.v = reshape(obj.calc.thres_v(:,j),ngrid,ngrid);
      else
        obj.print(['performing cloud segmentation for frame #' num2str(j)],1);
        r2b_midrange = obj.conf.seg.r2b_midrange;
        r2b = obj.curseg.r2b; r2b(~sm) = NaN;
        r2bmin = min(r2b_midrange(1),min(obj.curseg.r2b(obj.mfi.sm)));
        r2bmax = max(r2b_midrange(2),max(obj.curseg.r2b(obj.mfi.sm)));
        res = obj.conf.seg.r2b_resolution;
        resh = obj.conf.seg.r2b_hist_resolution;
        N = round((r2bmax-r2bmin)/res)+2;
        obj.curseg.thres.zz = r2bmin+(r2bmax-r2bmin-N*res)/2+res*(0:N);
        obj.curseg.thres.zzh = r2bmin-resh*1.5:resh:r2bmax+resh*2;
        NN = length(obj.curseg.thres.zzh);
        jmid_lo = find(obj.curseg.thres.zzh>=r2b_midrange(1),1);
        jmid_hi = find(obj.curseg.thres.zzh<=r2b_midrange(2),1,'last');
        obj.curseg.thres.X = tanh(bsxfun(@minus,obj.curseg.thres.zzh',obj.curseg.thres.zz)/res);
        obj.curseg.thres.ww0 = nan(N+1,ngrid,ngrid);
        obj.curseg.thres.ww = nan(N+1,ngrid,ngrid);
        obj.curseg.thres.bb = nan(ngrid);
        obj.curseg.thres.j0 = nan(ngrid);
        obj.curseg.thres.isouter = obj.thres.isouter;
        obj.curseg.thres.histn = zeros(NN,ngrid,ngrid);
        obj.curseg.thres.v0_mask = zeros(size(obj.curseg.sm));
        obj.curseg.thres.vlb = nan(ngrid); %presently not used
        obj.curseg.thres.vub = nan(ngrid); %presently not used
%         r2b1=r2b((obj.mfi.YY-sunp(1)).^2+(obj.mfi.XX-sunp(2)).^2<=15^2);
%         v_sun = mean(r2b1(~isnan(r2b1)));
        [~,sunpy_mfi]=min(abs(obj.mfi.yy-sunp(1)));
        [~,sunpx_mfi]=min(abs(obj.mfi.xx-sunp(2)));
        sunp_mfi = [sunpy_mfi;sunpx_mfi];
        allnn = zeros(NN,1);
        for ix=1:ngrid
          for iy=1:ngrid
            if sum(([obj.thres.yyg(iy+1);obj.thres.xxg(ix+1)]-sunp_mfi).^2) <= ...
                obj.conf.seg.rsun_mfi4outer^2
              obj.curseg.thres.isouter(iy,ix) = 0;
            end
%             if obj.curseg.clear_sun_flag==1 &&...
%                 abs(sunp_mfi(1)-obj.thres.yyg(iy+1))==min(abs(sunp_mfi(1)-obj.thres.yyg)) && ...
%                 abs(sunp_mfi(2)-obj.thres.xxg(ix+1))==min(abs(sunp_mfi(2)-obj.thres.xxg))
%               obj.curseg.thres.vub(iy,ix) = v_sun;
%             end
            z = r2b(obj.thres.yyg(iy):obj.thres.yyg(iy+2),obj.thres.xxg(ix):obj.thres.xxg(ix+2));
            z = z(:);
            jjy = obj.thres.yygb(iy)+1:obj.thres.yygb(iy+1);
            jjx = obj.thres.xxgb(ix)+1:obj.thres.xxgb(ix+1);
            if mean(isnan(z))>obj.conf.seg.tile_max_nan
              obj.curseg.thres.v0_mask(jjy,jjx) = -2;
            else
              z = z(~isnan(z));
              nn=hist(z,obj.curseg.thres.zzh)';
              w = wlsq(obj.curseg.thres.X,nn,[],...%.1+sqrt(mean(nn)./max(nn,1)),...
                obj.conf.seg.alpha*ones(1,N+1),[],[],...
                -obj.curseg.thres.X,nn*0,ones(1,N+1),0);
              yy = obj.curseg.thres.X*w;
              scale = sum(nn)/prod(2*obj.mfi.sz/ngrid)/sum(yy);
              yy = yy*scale; w = w*scale; nn = nn*scale;
              j0 = otsu_local_min(yy);
              if min(sum(yy(1:j0-1)),sum(yy(j0+1:end)))<sum(yy)*obj.conf.seg.outside_mass
                j0 = NaN;
              end
              obj.curseg.thres.j0(iy,ix) = j0;
              obj.curseg.thres.ww0(:,iy,ix) = w;
              obj.curseg.thres.histn(:,iy,ix) = nn;
              allnn = allnn+nn;
              if isnan(j0)
                obj.curseg.thres.v0_mask(jjy,jjx) = -1;
              else
                obj.curseg.thres.v0_mask(jjy,jjx) = obj.curseg.r2b(jjy,jjx)>=obj.curseg.thres.zzh(j0);
              end
              %repeat the fitting, but with high penalties at both ends
              b = max(nn);
              nnavg = mean(nn);
              s1 = sum(nn)*obj.conf.seg.outside_mass;
              j1 = max(1,min(jmid_lo,find(cumsum(nn)>=s1,1))-obj.conf.seg.n_relax_outside);
              j2 = max(1,min(NN+1-jmid_hi,find(cumsum(nn(end:-1:1))>=s1,1))-obj.conf.seg.n_relax_outside);
              nn(1:j1) = max(nn(1:j1),b+(nnavg-b)*(0:j1-1)'/j1);
              nn(NN:-1:NN+1-j2) = max(nn(NN:-1:NN+1-j2),b+(nnavg-b)*(0:j2-1)'/j2);
              w = -wlsq(obj.curseg.thres.X,b-nn,[],...
                obj.conf.seg.alpha*ones(1,N+1),[],[],...
                [-obj.curseg.thres.X;obj.curseg.thres.X],[nn*0;nn*0+b],ones(1,N+1),0);
              obj.curseg.thres.ww(:,iy,ix) = w;
              obj.curseg.thres.bb(iy,ix) = b;
            end
          end
        end
        
        if ~any(~isnan(obj.curseg.thres.j0(:)))
          %uniform sky, decide if uniformly cloudy or clear
          %allnn = allnn/sum(allnn);
          %if sum(obj.curseg.thres.zzh'.*allnn)>mean(obj.conf.seg.r2b_midrange)
          if obj.curseg.clear_sun_flag<=2
            j0 = 1; %max(1,find(cumsum(allnn)<=obj.conf.seg.extreme_outside_mass,1,'last')-1);
          else
            %disp(['WARNING: uniformly clear sky in frame #' num2str(j) ', this is unlikely']);
            j0 = NN; %min(NN,find(cumsum(allnn)>=1-obj.conf.seg.extreme_outside_mass,1)+1);
          end
          obj.curseg.thres.j0 = zeros(obj.conf.seg.ngrid)+j0;
          obj.curseg.thres.j0(~obj.thres.have_v) = NaN;
          obj.curseg.thres.v = zeros(obj.conf.seg.ngrid)+obj.curseg.thres.zzh(j0);
          obj.curseg.thres.v(~obj.thres.have_v) = NaN;
        else
          obj.curseg.thres.v = vmlSmoothThresMap(obj.curseg.thres,obj.thres.have_v,obj.conf.seg);
        end
        
        obj.calc.thres_v(:,j) = obj.curseg.thres.v(:);
      end
      
      [YYg,XXg] = ndgrid(obj.thres.yyg,obj.thres.xxg);
      VVg = XXg;
      VVg(2:end-1,2:end-1) = obj.curseg.thres.v;
      VVg(1,:) = VVg(2,:); VVg(end,:) = VVg(end-1,:);
      VVg(:,1) = VVg(:,2); VVg(:,end) = VVg(:,end-1);
      [YY,XX] = ndgrid(1:length(obj.mfi.yy),1:length(obj.mfi.xx));
      obj.curseg.thres.vmfi = interp2(XXg,YYg,VVg,XX,YY);
      obj.curseg.thres.vmfi(~obj.mfi.sm) = NaN;
      
      obj.curseg.d2sun = (obj.mfi.YY-sunp(1)).^2+(obj.mfi.XX-sunp(2)).^2;
      obj.curseg.d2sun_im = double(255*obj.curseg.d2sun/norm(obj.curseg.d2sun));
      obj.curseg.cc = double(obj.curseg.r2b>obj.curseg.thres.vmfi);
      obj.curseg.cc(isnan(obj.curseg.r2b) | ~obj.curseg.sm) = NaN;
      if obj.curseg.clear_sun_flag<0 || obj.curseg.clear_sun_flag>=3
        obj.curseg.cc(obj.curseg.r2b<obj.curseg.thres.vmfi & ...
          obj.curseg.d2sun>=obj.conf.rsun_px_close^2 & obj.curseg.d2sun<=obj.conf.rsun_px^2) = 0;
      end
      
      % amrollah: reflection
        obj.curseg.reflect= obj.curseg.cc .* double(rgb2gray(obj.curseg.x)-150);
        obj.curseg.reflect_fact = sum(sum(obj.curseg.reflect));
%         figure;imshow(reflect);
      
      obj.curseg.j = j;      
    end
    
% end cloud segmentation

%% process sequence

    function [jstart,jend] = find_daylight_range(obj,extend)
      if nargin<2, extend=0; end
      jj = obj.data.IrrClear(:,2)>=max(obj.data.IrrClear(:,2))*obj.conf.daylight_thres;
      jstart = obj.getidx(obj.data.IrrClear(find(jj,1),1)-extend*max(obj.conf.pred.tpred)/86400);
      jend = obj.getidx(obj.data.IrrClear(find(jj,1,'last'),1));
    end

    function xof = get_xof(obj)
        % for motion vector
      xof = double(obj.curseg.x(:,:,1));
      xof = xof-mean(xof(obj.curseg.sm));
      xof = xof/std(xof(obj.curseg.sm))+1j*obj.curseg.cc;
      xof(~obj.curseg.sm) = NaN;
    end
    
    function update_scal(obj,fn,upd,vupd)
        % update weigh of predict adaptation
      wold = obj.curpred.(['w' fn]);
      v = obj.curpred.(fn);
      wnew = repmat(obj.conf.scal.w_newest(:),size(wold,1)/length(obj.conf.scal.w_newest),2);
      wold = wold.*(1-wnew);
      wold(~upd) = max(wold(~upd),min(1,wnew(~upd)./(1-wnew(~upd))));
      wnew(~upd) = 0;
      v(upd) = (wold(upd).*v(upd)+wnew(upd).*vupd)./(wold(upd)+wnew(upd));
      obj.curpred.(fn) = v;
      obj.curpred.(['w' fn]) = wold+wnew;
    end
    
    function calc_cur_avgcc(obj)
        % calc cloud coverage of plant at frame j
      j = obj.curmot.j;
      assert(j==obj.curseg.j,'can do calc_cur_avgcc only if curmot and curmot point to the same frame');
      obj.calc.avgcc_sun(j) = obj.calc_avgcc_sun;
      obj.calc.avgcc(:,j) = obj.calc_avgcc_plant;
      
      %adapt scaling avgcc to irr / power
      irr = obj.getIrr(j)/obj.getIrrClear(j);
      upd = repmat([(obj.calc.avgcc_sun(j)<=obj.conf.scal.thres && obj.curseg.clear_sun_flag==4) ...
        (obj.calc.avgcc_sun(j)>=1-obj.conf.scal.thres && obj.curseg.clear_sun_flag==1)],...
        length(obj.conf.scal.w_newest),1);
      obj.update_scal('sun_avgcc2irr',upd,irr);

      p = obj.getP(j)/obj.getIrrClear(j);
      upd = repmat([obj.calc.avgcc(:,j)<=obj.conf.scal.thres obj.calc.avgcc(:,j)>=1-obj.conf.scal.thres],...
        [1 1 length(obj.conf.scal.w_newest)]);
      obj.update_scal('plant_avgcc2p',reshape(permute(upd,[3 1 2]),numel(upd)/2,2),p);

      obj.calc.sun_avgcc2irr(:,j) = obj.curpred.sun_avgcc2irr(:);
      obj.calc.plant_avgcc2p(:,j) = obj.curpred.plant_avgcc2p(:);
    end
    
    function update_pred(obj)
      j = obj.curmot.j;
      cc0 = obj.curmot.cc; cc0(isinf(obj.curmot.age_cc)) = nan;
      [j_p1,j_p2] = vmlMakeIndexPair(length(obj.conf.plantproj.extensions)*obj.conf.pred.n_cc2pow,...
        length(obj.conf.scal.w_newest),length(obj.conf.plantproj.heights));
      [j_irr1,j_irr2] = vmlMakeIndexPair(length(obj.conf.pred.rsun_px),length(obj.conf.scal.w_newest));
      irc0 = obj.getIrrClear(j);
      for i_t=1:length(obj.conf.pred.tpred)
        tpred = obj.conf.pred.tpred(i_t);
        irc1 = obj.getIrrClear(obj.data.ti(j)+tpred/86400);
        obj.curpred.p_persist = obj.getP(j)*irc1/irc0;
        obj.curpred.irr_persist = obj.getIrr(j)*irc1/irc0;
        p_true = obj.getP(obj.data.ti(j)+tpred/86400);
        irr_true = obj.getIrr(obj.data.ti(j)+tpred/86400);
        
        if tpred<=obj.conf.pred.tpred_persist
          obj.curpred.p_agg(i_t) = obj.curpred.p_persist;
          obj.curpred.irr_agg(i_t) = obj.curpred.irr_persist;
        else
          obj.curpred.p(i_t).v = nan+obj.curpred.p(i_t).sc;
          obj.curpred.p(i_t).v(end) = obj.curpred.p_persist;
          obj.curpred.irr(i_t).v = nan+obj.curpred.irr(i_t).sc;
          obj.curpred.irr(i_t).v(end) = obj.curpred.irr_persist;
          k_p = 0; k_irr = 0; bb=[];
          for i_v = 1:length(obj.conf.pred.n_mvec)
            n_mvec1 = obj.conf.pred.n_mvec(i_v);
            if any(isnan(obj.calc.mvec(:,j+1-n_mvec1))), break; end
            mvec1 = mean(obj.calc.mvec(:,j+1-n_mvec1:j),2);
            if norm(mvec1)<=obj.conf.mot.zero_tol
              p1 = obj.curpred.p_persist;
              irr1 = obj.curpred.irr_persist;
            else
              cc = vmlShift(obj.mfi,mvec1*tpred,cc0);
              [p1,bb] = obj.calc_avgcc_plant(obj.data.ti(j)+tpred/86400,cc,bb);
              p1 = vml_cc2pow(p1);
              p1 = irc1*(obj.curpred.plant_avgcc2p(j_p2,2).*(1-p1(j_p1))+obj.curpred.plant_avgcc2p(j_p2,1).*p1(j_p1));
              irr1 = 1-obj.calc_avgcc_sun(obj.data.ti(j)+tpred/86400,cc);
              irr1 = irc1*(obj.curpred.sun_avgcc2irr(j_irr2,2).*(1-irr1(j_irr1))+obj.curpred.sun_avgcc2irr(j_irr2,1).*irr1(j_irr1));
            end
            obj.curpred.p(i_t).v(k_p+(1:length(j_p1))) = p1;
            k_p = k_p+length(j_p1);
            obj.curpred.irr(i_t).v(k_irr+(1:length(j_irr1))) = irr1;
            k_irr = k_irr+length(j_irr1);
          end
          fupd = obj.conf.pred.upd_score;
          if mean(~isnan(obj.curpred.p(i_t).v))<obj.conf.pred.min_nexperts
            obj.curpred.p_agg(i_t) = NaN;
          else
            jj = ~isnan(obj.curpred.p(i_t).v);
            sc1 = obj.curpred.p(i_t).sc(jj); sc1=sc1-min(sc1);
            w = exp(-sc1/quantile(sc1,obj.conf.pred.quantil_score)); w = w/sum(w);
            obj.curpred.p_agg(i_t) = dot(w,obj.curpred.p(i_t).v(jj));
            err1 = abs(obj.curpred.p(i_t).v(jj)-p_true);
            obj.curpred.p(i_t).sc(jj) = (1-fupd)*obj.curpred.p(i_t).sc(jj)+fupd*err1;
            obj.curpred.p(i_t).sc(~jj) = (1-fupd)*obj.curpred.p(i_t).sc(~jj)+fupd*mean(err1);
          end
          obj.calc.p_pred(i_t,j) = obj.curpred.p_agg(i_t);
          
          if mean(~isnan(obj.curpred.irr(i_t).v))<obj.conf.pred.min_nexperts
            obj.curpred.irr_agg(i_t) = NaN;
          else
            jj = ~isnan(obj.curpred.irr(i_t).v);
            sc1 = obj.curpred.irr(i_t).sc(jj); sc1=sc1-min(sc1);
            w = exp(-sc1/quantile(sc1,obj.conf.pred.quantil_score)); w = w/sum(w);
            obj.curpred.irr_agg(i_t) = dot(w,obj.curpred.irr(i_t).v(jj));
            err1 = abs(obj.curpred.irr(i_t).v(jj)-irr_true);
            obj.curpred.irr(i_t).sc(jj) = (1-fupd)*obj.curpred.irr(i_t).sc(jj)+fupd*err1;
            obj.curpred.irr(i_t).sc(~jj) = (1-fupd)*obj.curpred.irr(i_t).sc(~jj)+fupd*mean(err1);
          end
          obj.calc.irr_pred(i_t,j) = obj.curpred.irr_agg(i_t);
        end
        
        obj.calc.p_abserr(i_t,j) = abs(obj.curpred.p_agg(i_t)-p_true);
        obj.calc.p_persist_abserr(i_t,j) = abs(obj.curpred.p_persist-p_true);
        obj.calc.irr_abserr(i_t,j) = abs(obj.curpred.irr_agg(i_t)-irr_true);
        obj.calc.irr_persist_abserr(i_t,j) = abs(obj.curpred.irr_persist-irr_true);
      end
    end
    
    function startframe(obj, j)
      obj.print(['starting at frame #' num2str(j)],1);
      obj.loadframe(j);
      obj.curmot.x = obj.get_xof;
      obj.curmot.j = j;
      obj.curmot.v = [0 0];
      obj.curmot.count = 1;
      obj.curmot.cc = obj.curseg.cc;
      obj.curmot.cc(isnan(obj.curmot.cc)) = 0;
      obj.curmot.w_cc = ~isnan(obj.curseg.cc)*obj.conf.mot.w_cc_newest;
      obj.curmot.age_cc = ones(size(obj.curseg.sm));
      obj.curmot.age_cc(isnan(obj.curseg.cc)) = inf;
      obj.calc_cur_avgcc;
    end
    
    function nextframe(obj)
      assert(~isempty(obj.curmot.x),'nextframe called without startframe');
      obj.print(['processing frame #' num2str(obj.curmot.j+1)],1);
      curmot1 = obj.curmot;
      j = curmot1.j+1;
      curmot1.j = j;
      curmot1.count = curmot1.count+1;
      obj.loadframe(j);
      curmot1.x = cat(3,curmot1.x,obj.get_xof);
      if size(curmot1.x,3)>obj.conf.mot.nframes
        curmot1.x = curmot1.x(:,:,end-obj.conf.mot.nframes+1:end);
      end
      tt = obj.data.ti(j-size(curmot1.x,3)+1:j);
      tt = (tt-tt(end))*86400;
      while any(abs(tt-(tt(1)+mean(diff(tt))*(0:length(tt)-1)))>1)
        curmot1.x = curmot1.x(:,:,2:end);
        tt = tt(2:end);
      end
      if all(~isnan(obj.calc.mvec(:,j))) 
        curmot1.v = obj.calc.mvec(:,j);
      else
        curmot1.v = obj.optical_flow(curmot1);
        obj.calc.mvec(:,j) = curmot1.v(:);
      end
      [curmot1.cc,curmot1.w_cc,curmot1.age_cc] = vmlShift(obj.mfi,...
        curmot1.v*(tt(end)-tt(end-1)),curmot1.cc,curmot1.w_cc,curmot1.age_cc);
      sm1 = ~isnan(obj.curseg.cc);
      wold_discounted = curmot1.w_cc(sm1)*(1-obj.conf.mot.w_cc_newest);
      curmot1.w_cc(sm1) = wold_discounted+obj.conf.mot.w_cc_newest;
      curmot1.cc(sm1) = (obj.curseg.cc(sm1)*obj.conf.mot.w_cc_newest+...
        curmot1.cc(sm1).*wold_discounted)./curmot1.w_cc(sm1);
      age_cc1 = ones(size(obj.curseg.sm));
      age_cc1(isnan(obj.curseg.cc)) = inf;
      curmot1.age_cc = min(curmot1.age_cc+1,age_cc1);
            
      obj.curmot = curmot1;
      obj.calc_cur_avgcc;
      obj.update_pred;
    end
    
    function process(obj,oldpred,jstart,jend)
      if nargin<2, oldpred=[]; end
      if nargin<4
        [jstart1,jend] = obj.find_daylight_range(1);
        if nargin<3, jstart = jstart1; end
      end
      obj.curpred = [];
      n_sun_avgcc2irr = length(obj.conf.scal.w_newest);
      n_plant_avgcc2p = length(obj.conf.plantproj.heights)*length(obj.conf.scal.w_newest);
      if isempty(oldpred)
        obj.curpred.sun_avgcc2irr = zeros(n_sun_avgcc2irr,2);
        obj.curpred.wsun_avgcc2irr = zeros(n_sun_avgcc2irr,2);
        obj.curpred.plant_avgcc2p = zeros(n_plant_avgcc2p,2);
        obj.curpred.wplant_avgcc2p = zeros(n_plant_avgcc2p,2);
      else
        assert(size(oldpred.sun_avgcc2irr,1)==n_sun_avgcc2irr,'size of oldpred.sun_avgcc2irr does not match');
        assert(size(oldpred.plant_avgcc2p,1)==n_plant_avgcc2p,'size of oldpred.plant_avgcc2p does not match');
        w0 = min(1,obj.conf.scal.w_newest(:)./(1-obj.conf.scal.w_newest(:)));
        obj.curpred.sun_avgcc2irr = oldpred.sun_avgcc2irr;
        obj.curpred.wsun_avgcc2irr = repmat(w0,1,2);
        obj.curpred.plant_avgcc2p = oldpred.plant_avgcc2p;
        obj.curpred.wplant_avgcc2p = repmat(w0,length(obj.conf.plantproj.heights),2);
      end
      nexperts = n_plant_avgcc2p*length(obj.conf.plantproj.extensions)*...
        obj.conf.pred.n_cc2pow*length(obj.conf.pred.n_mvec)+1;
      nexperts_sun = n_sun_avgcc2irr*length(obj.conf.pred.rsun_px)*...
        length(obj.conf.pred.n_mvec)+1;
      for i=1:length(obj.conf.pred.tpred)
        obj.curpred.p(i).sc = zeros(nexperts,1);
        obj.curpred.irr(i).sc = zeros(nexperts_sun,1);
      end
      obj.startframe(jstart);
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
      rsun_px = obj.conf.rsun_px; ls = '-';
      if isfield(obj.calc,'clear_sun_flag') && length(obj.calc.clear_sun_flag)>=j && ...
          isnan(obj.calc.clear_sun_flag(j))
        if obj.calc.clear_sun_flag(j)==1, ls = '--';
        elseif obj.calc.clear_sun_flag(j)==2, rsun_px = obj.conf.rsun_px_close;
        end
      end
      xx = yx(2)+cos(tt)*rsun_px;
      xx(xx<1 | xx>obj.oi.sz(2)) = NaN;
      yy = yx(1)+sin(tt)*rsun_px;
      yy(yy<1 | yy>obj.oi.sz(1)) = NaN;
      if ~isempty(marker)
        line(yx(2),yx(1),'color',col,'marker',marker); 
      end
      if ~isempty(marker2) && isfield(obj.calc,'yx_sun_detected') && ...
          size(obj.calc.yx_sun_detected,2)>=j && all(~isnan(obj.calc.yx_sun_detected(:,j))) &&...
          all(obj.calc.yx_sun_detected(:,j)>0)
        line(obj.calc.yx_sun_detected(2,j),obj.calc.yx_sun_detected(1,j),'color',col,'marker',marker2); 
      end
      line(xx,yy,'color',col,'linestyle',ls);
    end
    
    function showframe(obj,j,show_sun_or_skymask,h_fig)
      %show the frame, part which does not belong to the sky is
      %greened out
      if nargin<2, j=obj.curseg.j; end
      if nargin<3, show_sun_or_skymask = 0; end
      if nargin<4, h_fig = figure; else h_fig=clf(h_fig,'reset'); end
      if ~isscalar(show_sun_or_skymask), sm = show_sun_or_skymask;
      elseif show_sun_or_skymask>=2, sm = obj.skymask_wo_sun(j);
      else sm = obj.oi.sm; 
      end
      x = vmlColorify(obj.imread(j),~sm,2,64);
      image(x); axis off;
      set(h_fig,'KeyPressFcn',@(h_obj,evt) changeframe(evt));
      if isscalar(show_sun_or_skymask) && show_sun_or_skymask==1
        obj.plotSun(j,'g','.','x');
      end
      
      title([datestr(obj.data.ti(j),'HH:MM:SS') ' (#' num2str(j) ')']);
      function changeframe(evt)
          if numel(evt.Modifier)==1 && strcmp('shift',evt.Modifier{1})
              stepSize = 10;
          else
              stepSize = 1;
          end
          if strcmp(evt.Key, 'rightarrow')
              showframe(obj,min(j+stepSize,size(obj.data.ti,2)),show_sun_or_skymask,h_fig)
          elseif strcmp(evt.Key, 'leftarrow')
              showframe(obj,max(1,j-stepSize),show_sun_or_skymask,h_fig)
          end
      end
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
      title(obj.uidata.strtitle,'interpreter','none');
      datetick;
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
      datetick;
    end
    
    function plotirr(obj)
      %plot the power profiles
      plot(obj.data.ti,obj.getIrrClear(obj.data.ti),'r',...
        obj.data.ti,obj.getIrr(obj.data.ti),'b.');
      ylim = get(gca,'ylim'); ylim(1) = 0;
      set(gca,'ylim',ylim);
      grid on;
      datetick;
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
      title(obj.uidata.strtitle,'interpreter','none');
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
      title(obj.uidata.strtitle,'interpreter','none');
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
    
    function ShowSeg(obj,j)
      if nargin<2, 
        j = obj.curseg.j;
        assert(~isempty(j),'j not specified (previous computation was possibly aborted)');
      else
        obj.loadframe(j);
      end
      figno = 1322; obj.newfig(figno,['Cloud thresholds, frame #' num2str(j) ' / ' datestr(obj.data.ti(j))]);
      cb_text = {'r2b','glare','thres0','thres','thres val','contour'};
      hcb = zeros(1,length(cb_text));
      htxt = uicontrol('Style','text','Position',[10 5 650 16]);
      for i = 1:length(hcb)
        hcb(i) = uicontrol('Style','checkbox','Value',i==length(hcb),...
          'Position',[i*85 24 80 16],'String',cb_text{i});
      end
      auxplot_shown = [-1 -1];
      set(gcf,'toolbar','figure');
      set(gcf,'KeyPressFcn',@(h_obj,evt) changeframe(evt));
      set(hcb,'Callback',{@showSeg_upd}); %,hcb,hax
      showSeg_upd;
      
      function changeframe(evt)
          if numel(evt.Modifier)==1 && strcmp('shift',evt.Modifier{1})
              stepSize = 10;
          else
              stepSize = 1;
          end
          if strcmp(evt.Key, 'rightarrow')
              ShowSeg(obj,min(j+stepSize,size(obj.data.ti,2)));
          elseif strcmp(evt.Key, 'leftarrow')
              ShowSeg(obj,max(1,j-stepSize));
          end
      end
      
      function showSeg_click(varargin)
        p = get(gca,'currentpoint');
        [~,jx] = min(abs(p(1,1)-obj.mfi.xx));
        [~,jy] = min(abs(p(1,2)-obj.mfi.yy));
        ix = max(1,min(obj.conf.seg.ngrid,1+floor((p(1,1)-1)/(obj.oi.sz(2)-1)*obj.conf.seg.ngrid)));
        iy = max(1,min(obj.conf.seg.ngrid,1+floor((p(1,2)-1)/(obj.oi.sz(1)-1)*obj.conf.seg.ngrid)));
        set(htxt,'String',['xy = [' num2str(obj.mfi.xx(jx)) ', ' num2str(obj.mfi.yy(jy)) ...
           '] RGB = (' strjoin(cellfun(@num2str,squeeze(num2cell(obj.curseg.x(jy,jx,:))),'UniformOutput',0)',', ') ...
           ') R2B = ' num2str(obj.curseg.r2b(jy,jx),'%.2f') ' / ' ...
           'xytile = [' num2str(ix) ', ' num2str(iy) ']']);
         if any([ix iy]~=auxplot_shown) && isfield(obj.curseg.thres,'X')
           nn = obj.curseg.thres.histn(:,iy,ix);
           zzh = obj.curseg.thres.zzh;
           yy0 = obj.curseg.thres.X*obj.curseg.thres.ww0(:,iy,ix);
           yy = obj.curseg.thres.X*obj.curseg.thres.ww(:,iy,ix)+obj.curseg.thres.bb(iy,ix);
           j0 = obj.curseg.thres.j0(iy,ix);
           figure(13281); clf; hold on;
           bar(zzh,nn,1); 
           plot(zzh,max(0,yy0),'r',zzh,max(0,yy),'g--','linewidth',2);
           if ~isnan(j0), plot(zzh(j0),max(0,yy0(j0)),'ro','linewidth',2); end
           plot(obj.curseg.thres.v(iy,ix),0,'gd','linewidth',2);
           hold off; grid on;
         end
      end
      
      function showSeg_upd(varargin)
        figure(figno);
        vv = get(hcb,'value');
        if ~vv{1}, x = obj.curseg.x;
        else
          x = obj.curseg.r2b; x = x-min(x(:));
          x = repmat(uint8(x/max(x(:))*192),[1 1 3]);
        end
        x = vmlColorify(x,~obj.mfi.sm,1:3,-255);
        if vv{2}
            x = vmlColorify(x,obj.curseg.glare,2,80); 
%           x = vmlColorify(x,isnan(obj.curseg.r2b) & ~isnan(obj.curseg.sm),2,80); 
        end
        if vv{4}
          x = vmlColorify(x,obj.curseg.cc==0,1,-40);
          x = vmlColorify(x,obj.curseg.cc==1,1,40);
        elseif vv{3} && isfield(obj.curseg.thres,'v0_mask')
          x = vmlColorify(x,obj.curseg.thres.v0_mask==-2,1:3,-80);
          x = vmlColorify(x,obj.curseg.thres.v0_mask==-1,1:3,40);
          x = vmlColorify(x,obj.curseg.thres.v0_mask==0,1,-40);
          x = vmlColorify(x,obj.curseg.thres.v0_mask==1,1,40);
        end
        ih = image(obj.mfi.xx,obj.mfi.yy,x); axis off;
        title([datestr(obj.data.ti(j),'HH:MM:SS') ' (#' num2str(j) ')']);
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
        if vv{5} && length(unique(obj.curseg.thres.vmfi(~isnan(obj.curseg.thres.vmfi))))>1
          hold on;
          [~,ih] = contour(obj.mfi.XX,obj.mfi.YY,obj.curseg.thres.vmfi,20);
          set(gca,'clim',[min(obj.curseg.thres.vmfi(:)) max(obj.curseg.thres.vmfi(:))]);
          hold off;
          set(ih,'ButtonDownFcn',@showSeg_click)
          pos = get(gca,'Position');
          colorbar('vert');
          set(gca,'Position',pos);
        elseif vv{6} && length(unique(obj.curseg.cc(~isnan(obj.curseg.cc))))>1
          hold on;
          [~,ih] = contour(obj.mfi.XX,obj.mfi.YY,obj.curseg.cc,[.5 .5],'g','linewidth',1);
          hold off;
          set(ih,'ButtonDownFcn',@showSeg_click)
        end
        set(gca,'xlim',[0.5 obj.oi.sky_area(4)-obj.oi.sky_area(3)+1.5]);
        set(gca,'ylim',[0.5 obj.oi.sky_area(2)-obj.oi.sky_area(1)+1.5]);
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
    
    function yx = sunpos_im(obj,t)
      %position of the sun in image pixels
      yx = obj.camworld2im(obj.calib.Rext'*obj.sunpos_realworld(t));
    end
    
    function yx = sunpos_mfi(obj,t)
      %position of the sun in mfi image pixels
      yx = round(obj.sunpos_im(t).*(obj.mfi.sz./obj.oi.sz)');
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
      if nargin<3
        cc = obj.curmot.cc; 
        cc(isinf(obj.curmot.age_cc)) = nan;
        y = mean(cc(~isnan(cc) & obj.curseg.d2sun<=obj.conf.pred.rsun_px(3)^2)); %amri added the first ~isnan
      else
        sunp = obj.sunpos_im(t);
        d2sun = (obj.mfi.YY-sunp(1)).^2+(obj.mfi.XX-sunp(2)).^2;
        y = zeros(length(obj.conf.pred.rsun_px),1);
        for i=1:length(y)
          y(i) = mean(cc(~isnan(cc) & d2sun<=obj.conf.pred.rsun_px(i)^2));%amri added the first ~isnan
        end
      end
    end
    
    function [y, bb] = calc_avgcc_plant(obj,t,cc,bb)
      if nargin<3, 
        t = obj.curmot.j; 
        cc = obj.curmot.cc; 
        cc(isinf(obj.curmot.age_cc)) = nan;
        ee = 0;
      else
        ee = obj.conf.plantproj.extensions;
      end
      has_bb = nargin>=4 && ~isempty(bb);
      hh = obj.conf.plantproj.heights;
      y = zeros(length(ee)*length(hh),1);
      k = 0;
      j1 = ~isnan(obj.mfi.ppx(:,1)) & ~isnan(cc(:));% amri added the secobd part
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
    
    function sm = skymask_wo_sun(obj,j,iminfo,rsun_px)
      %sky mask without sun, frame j (taking out a circle around 
      %the sun position)
      if nargin<3, iminfo = obj.oi; end
      if nargin<4, rsun_px = obj.conf.rsun_px; end
      yx = sunpos_im(obj,j);
      sm = iminfo.sm;
      sm((iminfo.XX-yx(2)).^2+(iminfo.YY-yx(1)).^2<=rsun_px^2) = false;
    end
    
    function w = weights4oflow(obj,j,mvec1,iminfo)
      %compute sky mask for the part of the sky where the clouds
      %come from
      if nargin<2, j=obj.curmot.j; end
      if nargin<3, mvec1 = obj.curmot.v; end
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
      obj.print(['#frames = ' num2str(length(frame_nos)) ', pixel RMSE = ' ...
        num2str(sqrt(mean(sum((pp-obj.calib.kabsch.im_model).^2,1))),'%.1f') ...
        ', RMSE on 1 m unit sphere = ' ...
        num2str(1000*sqrt(mean(sum((R*PP-QQ).^2,1))),'%.1f') ' mm'],1);
    end
    
% end sun detection

%% motion vectors / optical flow

    function v = optical_flow(obj, curmot1)
      if nargin<2, curmot1 = obj.curmot; end
      j = curmot1.j;
      obj.print(['computing motion vector for frame #' num2str(j)],1);
      if curmot1.count<=2
        curmot1.v = vmlOpticalFlow(obj.mfi,obj.conf.oflow,curmot1.x,...
          obj.data.ti(j-size(curmot1.x,3)+1:j),curmot1.v);
      end
      w = obj.weights4oflow;
      v = vmlOpticalFlow(obj.mfi,obj.conf.oflow,curmot1.x,...
        obj.data.ti(j-size(curmot1.x,3)+1:j),curmot1.v,w);
    end
  
% end motion vectors
  
  end 
end
