classdef vmlSeq < handle
  properties
    X;              %frames
    files;          %file names
    folder;         %folder
    dti;            %date - time for images
    dtm;            %date - time for measurements 
    P;              %PV plant power
    GHI;            %GHI measurement
    Temp;           %temperature
    conf;           %configuration struct
    sky_mask;       %sky mask
    sky_mask_border;%sky mask border
    sky_area;       %smallest rectangle that fully contains the sky
    sphere_px;      %pixels on the unit sphere
    plane_px;       %pixels on the plane at height 1
    calib;          %camera calibration
    strday;         %day as a string
    dplane;         %min pixel distance on plane
    mvec;           %motion vectors on plane
    printlevel = 0; %0=silent, 1=summarizing, 2=info, 3=verbose
  end
  
  methods
    
    function obj = vmlSeq(month_day,hours,conf)
      if nargin<2, hours = []; end
      if nargin<3, conf = evalin('base','VMLCONF'); end
      obj.conf = conf;
      [obj.files, obj.folder, obj.dti] = vmlListImgFiles(month_day, hours, conf);
      
      x1=load(conf.sky_mask);
      fn = fieldnames(x1);
      obj.sky_mask = x1.(fn{1});
      obj.sky_mask = obj.sky_mask(:,:,1);
      
      obj.calib = [];
      x1=load(conf.calibration{1});
      fn = fieldnames(x1);
      obj.calib.intern = x1.(fn{1});
      x1=load(conf.calibration{2});
      fn = fieldnames(x1);
      obj.calib.Rext = x1.(fn{1});
      x1=load(conf.calibration{3});
      fn = fieldnames(x1);
      obj.calib.model3D = x1.(fn{1});
      
      height = obj.calib.intern.height;
      width = obj.calib.intern.width;
      
      [yy,xx] = ndgrid(1:height,1:width);
      obj.sphere_px = cam2world([yy(:)';xx(:)'],obj.calib.intern);
      obj.sky_mask(obj.sphere_px(3,:)>-sin(conf.minangle)) = false;
      
      obj.sphere_px = reshape(obj.sphere_px',[height width 3]);
      
      obj.sky_area = [max(0,find(any(obj.sky_mask,2),1)-1) ...
        min(height,find(any(obj.sky_mask,2),1,'last')+1) ...
        max(0,find(any(obj.sky_mask,1),1)-1) ...
        min(width,find(any(obj.sky_mask,1),1,'last')+1)];
      
      obj.sky_mask = obj.sky_mask(obj.sky_area(1):obj.sky_area(2),obj.sky_area(3):obj.sky_area(4));
      
      obj.sky_mask_border = false(size(obj.sky_mask));
      obj.sky_mask_border(1:end-1,:) = obj.sky_mask_border(1:end-1,:) | obj.sky_mask(2:end,:);
      obj.sky_mask_border(2:end,:) = obj.sky_mask_border(2:end,:) | obj.sky_mask(1:end-1,:);
      obj.sky_mask_border(:,1:end-1) = obj.sky_mask_border(:,1:end-1) | obj.sky_mask(:,2:end);
      obj.sky_mask_border(:,2:end) = obj.sky_mask_border(:,2:end) | obj.sky_mask(:,1:end-1);
      obj.sky_mask_border(obj.sky_mask) = false;
      
      obj.sphere_px = obj.sphere_px(obj.sky_area(1):obj.sky_area(2),obj.sky_area(3):obj.sky_area(4),:);
      obj.sphere_px = reshape(obj.sphere_px,numel(obj.sky_mask),3);
      obj.plane_px = -obj.sphere_px(:,1:2)./repmat(obj.sphere_px(:,3),[1 2]);
      obj.plane_px(~(obj.sky_mask | obj.sky_mask_border),:) = inf;
      
      plane_px1 = reshape(obj.plane_px,[size(obj.sky_mask) 2]);
      d1 = plane_px1(2:end,:,1)-plane_px1(1:end-1,:,1);
      d2 = plane_px1(:,2:end,2)-plane_px1(:,1:end-1,2);
      obj.dplane = min(min(abs(d1(~isnan(d1)))),min(abs(d2(~isnan(d2)))));
      assert(obj.dplane >= 1/max(height,width));
      
      x1=load(conf.measured_data);
      fn = fieldnames(x1);
      x1 = x1.(fn{1});
      j1 = max(1,find(x1.timestamp>=min(obj.dti),1)-1);
      j2 = min(length(x1.timestamp),find(x1.timestamp<=max(obj.dti),1,'last')+1);
      obj.dtm = x1.timestamp(j1:j2);
      obj.P = x1.power(j1:j2);
      obj.GHI = x1.GHI(j1:j2);
      obj.Temp = x1.temperature(j1:j2);
      obj.strday = [num2str(conf.year) '_' num2str(month_day(1),'%02d')  '_' num2str(month_day(2),'%02d')];
      
      obj.mvec = zeros(2,length(obj.dti))+NaN;
    end
    
    function p = sunpos_realworld(obj,j)
      t = [];
      [t.year,t.month,t.day,t.hour,t.min,t.sec]=datevec(obj.dti(j));
      t.UTC = obj.calib.model3D.UTC;
      sun = sun_position(t, obj.calib.model3D.Location);
      [x,y,z] = sph2cart(sun.azimuth*pi/180,(90-sun.zenith)*pi/180,1);
      p = [x y z]';
    end
    
%     function p = sunpos_img(obj,j)
%       t = [];
%       [t.year,t.month,t.day,t.hour,t.min,t.sec]=datevec(obj.dti(j));
%       t.UTC = obj.calib.model3D.UTC;
%       sun = sun_position(t, obj.calib.model3D.Location);
%       [x,y,z] = sph2cart((90-sun.azimuth)*pi/180,(90-sun.zenith)*pi/180,1); 
%       p = fliplr(world2cam(obj.calib.extern.Rinv*[x;y;-z],...
%         obj.calib.intern)')-obj.sky_area([3 1])+1;
%     end

    function p = sunpos_img(obj,j)
      p = obj.calib.Rext'*obj.sunpos_realworld(j);
      p(3) = -p(3);
      p = fliplr(world2cam(p,obj.calib.intern)')-obj.sky_area([3 1])+1;
    end

    function sm = skymask_wo_sun(obj,j)
      p = sunpos_img(obj,j);
      sm = obj.sky_mask;
      [yy,xx] = ndgrid(1:size(sm,1),1:size(sm,2));
      sm((xx-p(1)).^2+(yy-p(2)).^2<=obj.conf.rsun_px^2) = false;
    end
    
    function x = imread(obj,j)
      x = imread([obj.folder obj.files{j} obj.conf.img_ext]);
      x = x(obj.sky_area(1):obj.sky_area(2),obj.sky_area(3):obj.sky_area(4),:);
    end
    
    function showframe(obj,j,with_sun)
      if nargin<3, with_sun = 0; end
      x = obj.imread(j);
      x = reshape(x,[numel(obj.sky_mask) 3]);
      if with_sun, sm = obj.sky_mask;
      else sm = obj.skymask_wo_sun(j); end
      x(~sm,2) = min(255,x(~sm,2) + 128);
      x = reshape(x,[size(obj.sky_mask) 3]);
      image(x); axis off;
      title([datestr(obj.dti(j),'HH:MM:SS') ' (#' num2str(j) ')']);
    end
    
    function p = detect_saturated_sun(obj,j)
      x = sum(obj.imread(j),3)/(3*255);
      x(~(obj.sky_mask | obj.sky_mask_border)) = 0.5;
      x = max(abs(cv.Sobel(x,'XOrder',1,'YOrder',0)),...
              abs(cv.Sobel(x,'XOrder',0,'YOrder',1)));
      x(~obj.sky_mask) = 0;
      x=cv.GaussianBlur(x,'KSize',obj.conf.gaussblur4sundect_px+[0 0]);
      [~,jmax] = max(x(:));
      [p,q] = ind2sub(size(obj.sky_mask),jmax);
      %figure(55);imagesc(x);hold on;plot(q,p,'gx');hold off
      p = [q+obj.sky_area(3)-1 p+obj.sky_area(1)-1];
    end
    
    function R = kabsch(obj,frame_nos)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.dti); 
      end
      pp = zeros(2,length(frame_nos));
      PP = zeros(3,length(frame_nos)); QQ = PP;
      for j = 1:length(frame_nos)
        pp(:,j) = obj.detect_saturated_sun(j)';
        PP(:,j) = cam2world(flipud(pp(:,j)),obj.calib.intern);
        QQ(:,j) = obj.sunpos_realworld(j);
      end
      PP(3,:) = -PP(3,:);
      [U,~,V] = svd(PP*QQ');
      R = V*diag([1 1 sign(det(V*U'))])*U';
      obj.calib.kabsch.frame_nos = frame_nos;
      obj.calib.kabsch.dt = obj.dti(frame_nos);
      obj.calib.kabsch.world_detect = PP;
      obj.calib.kabsch.realword = QQ;
      obj.calib.kabsch.im_detect = pp;
      PP1 = R'*QQ;
      PP1(3,:) = -PP1(3,:);
      obj.calib.kabsch.im_model = flipud(world2cam(PP1,obj.calib.intern));
      obj.calib.kabsch.R = R;
    end
    
    function show_ecalib(obj)
      if ~isfield(obj.calib,'kabsch') || isempty(obj.calib.kabsch)
        error('no external calibration was executed on this vmlSeq object');
      end
      figure(1329); clf;
      subplot(3,2,5:6);
      plot(obj.calib.kabsch.dt,obj.calib.kabsch.im_model(1,:),'b',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_model(2,:),'r',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_detect(1,:),'bo',...
        obj.calib.kabsch.dt,obj.calib.kabsch.im_detect(2,:),'ro');
      title(obj.strday,'interpreter','none');
      datetickzoom;
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@show_ecalib_upd)
      function out = show_ecalib_upd( ~, event_obj)
        figure(1329);
        [~,j] = min(abs(event_obj.Position(1)-obj.calib.kabsch.dt));
        subplot(3,2,1:4); 
        obj.showframe(obj.calib.kabsch.frame_nos(j),1);
        hold on; 
        plot(obj.calib.kabsch.im_model(1,:)-obj.sky_area(3)+1,...
          obj.calib.kabsch.im_model(2,:)-obj.sky_area(1)+1,'r',...
          obj.calib.kabsch.im_detect(1,:)-obj.sky_area(3)+1,...
          obj.calib.kabsch.im_detect(2,:)-obj.sky_area(1)+1,'ro');
        hold off;
        drawnow;
        out = datestr(event_obj.Position(1),'HH:MM:SS');
      end
    end
    
    function show(obj)
      figure(1324); clf;
      subplot(2,3,4:6);
      plot(obj.dtm,obj.P,'b',obj.dti,interp1(obj.dtm,obj.P,obj.dti),'b.');
      hold on; hmarks = plot(obj.dti(1:3),interp1(obj.dtm,obj.P,obj.dti(1:3)),'ro'); hold off;
      set(hmarks, 'XData', [], 'YData', []);
      title(obj.strday,'interpreter','none');
      datetickzoom;
      dcmobj = datacursormode(gcf);
      set(dcmobj,'UpdateFcn',@show_upd)
      function out = show_upd( ~, event_obj)
        figure(1324);
        jj = zeros(1,3);
        for i=0:2
          [~,j] = min(abs(event_obj.Position(1)-i*obj.conf.show_dt/86400-obj.dti));
          subplot(2,3,3-i); 
          obj.showframe(j);
          jj(i+1) = j;
        end
        set(hmarks, 'XData', obj.dti(jj), 'YData', interp1(obj.dtm,obj.P,obj.dti(jj)));
        drawnow;
        out = datestr(event_obj.Position(1),'HH:MM:SS');
      end
    end
    
    function v = mvec_filtered(obj, frame_nos)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.dti); 
      end 
      v = zeros(2,numel(frame_nos)) + NaN;
      for i = 1:numel(frame_nos)
        j = frame_nos(i);
        if j<2 || j>length(obj.dti), continue; end
        n = min(obj.conf.nfilter_mvec,j-1);
        v(:,i) = mean(obj.mvec(:,j-n+1:j),2);
      end
    end
    
    function optical_flow(obj, frame_nos, do_replace)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 2:length(obj.dti); 
      end
      if nargin<3, do_replace = 0; end
      for frame_no = frame_nos(:)'
        if frame_no<2 || frame_no>length(obj.dti) || (~do_replace && all(~isnan(obj.mvec(:,frame_no))))
          continue; 
        end
        if obj.printlevel>=1
          disp(['computing motion vector for frame #' num2str(frame_no)]);
        end
        x0 = obj.imread(frame_no-1);
        x0 = double(x0(:,:,1)); %./max(80,sum(x00,3)); %red channel over sum
        x1 = obj.imread(frame_no);
        x1 = double(x1(:,:,1)); %./max(80,sum(x1,3)); %red channel over sum
        sm = obj.skymask_wo_sun(frame_no);
        jorig = find(sm);
        j = jorig;
        %figure(10);image(x00);
        %figure(11);imagesc(x0);
        %figure(12);imagesc(x1);
        %       jj = 629:829;
        %       smj=obj.sky_mask+0;
        %       smj(obj.sky_mask)=obj.jsky_mask;
        %       jjj=smj(jj,jj);
        %       j = jjj(:);
        %       jorig = j;
        
        fit = mean((x0(j)-x1(j)).^2);
        dx = [0 0];
        dx1 = obj.mvec(:,frame_no-1)';
        if all(~isnan(dx1))
          j1=closest_in_plane(obj.plane_px,jorig,j,dx1,size(sm),0);
          fit1 = sum(((x0(j1)-x1(jorig)).*sm(j1)).^2)/sum(sm(j1));
          if fit1<fit, dx = dx1;j = j1; fit = fit1; end
          dx1 = obj.mvec_filtered(frame_no)';
          if frame_no>3 && all(~isnan(dx1))
            j1=closest_in_plane(obj.plane_px,jorig,j,dx1,size(sm),0);
            fit1 = sum(((x0(j1)-x1(jorig)).*sm(j1)).^2)/sum(sm(j1));
            if fit1<fit, dx = dx1;j = j1; fit = fit1; end
          end
        end
        scale = 25;
        iprev = [];
        while 1
          j1best = [];
          for i=0:3
            if ~any(i==iprev)
              dx1 = dx + scale*obj.dplane*(-1)^rem(i,2)*([0 1]+floor(i/2)*[1 -1]);
              j1=closest_in_plane(obj.plane_px,jorig,j,dx1,size(sm),0+(scale==1));
              fit1 = sum(((x0(j1)-x1(jorig)).*sm(j1)).^2)/sum(sm(j1));
              if fit1<fit
                dx1best = dx1; fit = fit1; j1best = j1; iprev1 = bitxor(i,1);
              end
              if obj.printlevel>=3
                disp(['    ' num2str(round(dx1/obj.dplane)) ' ' num2str(fit1)]);
              end
            end
          end
          if ~isempty(j1best)
            dx = dx1best; j = j1best; iprev = iprev1;
            if obj.printlevel>=3
              disp(['*** ' num2str(round(dx/obj.dplane)) ' ' num2str(fit)]);
            end
          else
            if scale==1, break; end
            iprev = [];
            scale = round(scale/5);
          end
        end
        obj.mvec(:,frame_no) = dx(:);
      end
    end
    
  end 
end
