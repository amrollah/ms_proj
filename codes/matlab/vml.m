classdef vml < vmlSeq  
    properties
    end
    
    methods
        function obj = vml(plantid,folder,hours,orig)
            if nargin<4, orig = false; end
            if nargin<3, hours = []; end
            obj@vmlSeq(plantid,folder,hours);
            
            if ~orig
                %clear sky irradiation
                tt = (obj.data.tmin:6/86400:obj.data.tmax)';
                [ClearSkyGHI,ClearSkyDNI,ClearSkyDHI,~,~] = pvl_clearsky_ineichen(pvl_maketimestruct(tt, ...
                obj.calib.UTC),obj.calib.model3D.Location);
                obj.data.IrrClear = [tt, ClearSkyGHI, ClearSkyDNI, ClearSkyDHI];

                [ti_ClearSkyGHI, ti_ClearSkyDNI, ti_ClearSkyDHI, zenith, azimuth] = pvl_clearsky_ineichen(pvl_maketimestruct(obj.data.ti, ...
                obj.calib.UTC),obj.calib.model3D.Location);
                obj.data.ti_elevation = 90 - zenith;
                obj.data.zenith = zenith;
                obj.data.ti_azimuth = azimuth;
                obj.data.ti_IrrClear = [obj.data.ti', ti_ClearSkyGHI, ti_ClearSkyDNI, ti_ClearSkyDHI];
                
                obj.calc.Irr = NaN(length(obj.data.ti),4);
                obj.calc.Irr(:,1) = obj.data.ti(:);
            end
        end
        
        function loadframe(obj,j,redo_thres_calc)
            if nargin<3, redo_thres_calc=0; end
            loadframe@vmlSeq(obj,j,redo_thres_calc);
            
            if any(~isnan(obj.calc.thres_v(:,j)))
                % amri: diffuse irradiance calcualtion
                % obj.calc.Irr = [time, Diffuse, Direct, Global]
                % Diffuse = GHI-DNI*cosd(Z)*(1|0)
                diffuse = obj.data.Irr(j,2)-obj.data.ti_IrrClear(j,3)*cosd(obj.data.zenith(j))*obj.sunFlagToCoef(j);
                obj.calc.Irr(j,2) = diffuse;
%                 obj.print(1,['Diffuse Irr: ' num2str(diffuse)]);
            end
            % reflection
            cloud_shine = obj.curseg.cc .* double(rgb2gray(obj.curseg.x)-150);
            cloud_shine(cloud_shine<0)=0;
            obj.calc.seg.cloud_shine_fact(j) = nansum(nansum(cloud_shine));
            obj.calc.seg.clouds_fact(j) = floor(100*nansum(obj.curseg.cc(:))/sum(obj.mfi.sm(:)));
            % figure;imshow(obj.curseg.cloud_shine);
        end
        
         function coef = sunFlagToCoef(obj,j,mode)
             if nargin<3
                 mode=1;
             end
             if mode==1
                 clear_sun_flag = obj.calc.seg.clear_sun_flag(j);
             else
                 clear_sun_flag=j;
             end
              switch lower(clear_sun_flag)
                case 1 % no visible sun, DNI~0
                  coef = 0;
                case 4 % complete star shape sun, DNI~1
                  coef = 1;
                otherwise
                  coef = NaN;
              end
         end
         
         function p = sunpos_realworld(obj,t)
          %position of the sun in the world coordinates
          if t<length(obj.data.ti), t = obj.data.ti(t); end
          t1 = [];
          [t1.year,t1.month,t1.day,t1.hour,t1.min,t1.sec]=datevec(t);
          t1.UTC = obj.calib.model3D.UTC+isdst(datetime(t,'ConvertFrom','datenum','TimeZone',obj.conf.TimeZone));
          sun = sun_position(t1, obj.calib.model3D.Location);
          [x,y,z] = sph2cart(sun.azimuth*pi/180,(90-sun.zenith)*pi/180,1);
          p = [x y z]';
         end
        
         function I = get_image(obj, filename)
            I=imread(filename);
            if ~isempty(obj.oi.sky_area)
                I = I(obj.oi.sky_area(1):obj.oi.sky_area(2),obj.oi.sky_area(3):obj.oi.sky_area(4),:);
            end
         end
         
         function y = get45Irr(obj,t,tpred)
          %get the irradiation data for time(s) t
          if nargin<2, t = obj.data.ti;
          elseif any(t<obj.data.ti(1)-1), t=obj.data.ti(t); 
          end
          if nargin>=3, t = t+tpred/86400; end
          y = interp1(obj.data.Irr(:,1),obj.data.Irr(:,3),t);
        end

        function y = getDiffuseIrr(obj,t,tpred)
          %get the irradiation data for time(s) t
          if nargin<2
              y = obj.calc.Irr(:,2);
              return;
          elseif any(t<obj.data.ti(1)-1), t=obj.data.ti(t); 
          end
          if nargin>=3, t = t+tpred/86400; end
          [~, idt] = arrayfun(@(ti) min(abs(obj.calc.Irr(:,1)-ti)),t);
          y = obj.calc.Irr(idt,2);
        end


        function plot45irr(obj)
          %plot the irradiation profiles
          plot(obj.data.ti,obj.getIrrClear(obj.data.ti),'bl', obj.data.ti,obj.getIrr(obj.data.ti),'r',...
            obj.data.ti,obj.get45Irr(obj.data.ti),'b.-', obj.data.ti,obj.getDiffuseIrrClear(obj.data.ti),'g',...
        obj.data.ti,obj.getDiffuseIrr(obj.data.ti),'c.');
          ylim = get(gca,'ylim'); ylim(1) = 0;
          set(gca,'ylim',ylim);
          grid on;
          datetickzoom;
        end

        function [d,tilted_DNI] = get45diffuse(obj,t,az,elev)
            if nargin<4
              az = .7;
              elev = 33.5;
            end
            plate_cord = repmat([deg2rad(az),deg2rad(elev),1],[length(t),1]);
            [px,py,pz] = sph2cart(plate_cord(:,1),plate_cord(:,2),plate_cord(:,3));
            plate_cord = [px,py,pz];
            [sx,sy,sz] = sph2cart(deg2rad(obj.data.ti_azimuth(t)),deg2rad(obj.data.ti_elevation(t)),ones(size(t))');
            sun_cords = [sx,sy,sz];
    %         figure; plot3([px(1),0],[py(1),0],[pz(1),0]);
    %         hold on; scatter3(sx,sy,sz,'b');
            nrm=sqrt(sum(abs(cross(sun_cords,plate_cord,2)).^2,2));
            angles = atan2d(nrm, dot(sun_cords,plate_cord,2));
            tilted_DNI = obj.data.ti_IrrClear(t,3).*max(0,cosd(angles));
            d = obj.get45Irr(t)-tilted_DNI';%.*obj.sunFlagToCoef;
        end

        function plot45diffuse(obj,az,elev,plotdiff)
            if nargin<4, plotdiff=0; end
              %plot the diffuse irradiance derived from 45 degree sensor
              if nargin==1
                  az = .7;
                  elev = 33.5;
              end
              [df,DNI] = obj.get45diffuse((1:length(obj.data.ti)),az,elev);

              other_df= obj.getIrr()'-obj.data.ti_IrrClear(:,3).*cosd(obj.data.zenith(:));
              h = plot(obj.data.ti,obj.getIrr,obj.data.ti,obj.getDiffuseIrrClear(obj.data.ti),'r',...
                obj.data.ti,df,'b.',obj.data.ti,DNI,'c.-',...
                obj.data.ti,obj.get45Irr((1:length(obj.data.ti))),'g--',...
                obj.data.ti,other_df,'k.-');
            legend(h, 'HGI', 'diffudse clear(HDI)', 'tilted diffuse', 'tilted effective DNI', 'tilted irradiation', 'diffuse analytic');
            if plotdiff
                h = plot(obj.data.ti,df,'-k',obj.data.ti,other_df,'-r',obj.data.ti,df-other_df','b.');
            end
              ylim = get(gca,'ylim'); ylim(1) = 0;
              set(gca,'ylim',ylim);
              grid on;
              datetickzoom;
        end
        
        function y = getDiffuseIrrClear(obj,t,tpred)
          %get the clear sky diffuse irradiation for time(s) t
          if nargin<2, t = obj.data.ti;
          elseif any(t<obj.data.ti(1)-1), t=obj.data.ti(t); 
          end
          if nargin>=3, t = t+tpred/86400; end
          y = interp1(obj.data.IrrClear(:,1),obj.data.IrrClear(:,4),t);
        end
        
        function plotdiffuse(obj)
          %plot the diffuse irradiance
          plot(obj.data.ti,obj.getDiffuseIrrClear(obj.data.ti),'r',...
            obj.data.ti,obj.getDiffuseIrr(obj.data.ti),'b.-');
          ylim = get(gca,'ylim'); ylim(1) = 0;
          set(gca,'ylim',ylim);
          grid on;
          datetickzoom;
        end
        
        function ShowD(obj, method)
          if nargin<2, method = 'showtraj'; end
          figno = 11; obj.newfig(figno,'Diffuse vs Global irradiation');
          ax(1) = subplot(2,4,3:4); obj.plotdiffuse;
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
            if ~isfield(obj.calc, 'seg') || length(obj.calc.seg.clear_sun_flag)<j || obj.calc.seg.clear_sun_flag(j)==0
                obj.loadframe(j);
                subplot(2,4,3:4); obj.plotdiffuse;
            end
            drawnow;
            if nargin<2, out=[]; 
            else out = {datestr(event_obj.Position(1),'HH:MM:SS'), ['clouds: ',num2str(obj.calc.seg.clouds_fact(j)),'%'], ['cloud-shine: ', num2str(obj.calc.seg.cloud_shine_fact(j))]};
            end
          end
        end
    end
    
end

