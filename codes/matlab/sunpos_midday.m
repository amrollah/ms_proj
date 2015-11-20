function [pos, zenith] = sunpos_midday(obj)
      %position of the sun in the world coordinates at midday
      hour = 11;
      min = 20;
      t = obj.gen_time_struct(1,[hour,min,0]);
      sun = sun_position(t, obj.calib.model3D.Location);
      
      last_zenith = 180;
      while sun.zenith<last_zenith         
          last_zenith = sun.zenith;  
          min = min + 1;
          if min==60
              hour = hour + 1;
              min = 0;
          end
          t = obj.gen_time_struct(1,[hour,min,0]);
          sun = sun_position(t, obj.calib.model3D.Location);
      end
      [x,y,z] = sph2cart(sun.azimuth*pi/180,(90-sun.zenith)*pi/180,1);
      p = [x y z]';
      pos = obj.camworld2im(obj.ext_calib.R'*p);
      zenith = last_zenith;
      %fprintf('midday time: %d:%d:00   min_zenith=%f\n', hour, min, last_zenith);
end