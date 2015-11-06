function pos = sunpos_midday(obj)
    load([obj.conf.datafolder 'calib.mat']);
      %position of the sun in the world coordinates at 12:00
      date = textscan(obj.strtitle,'%f','delimiter','_');
      t.year = date{1,1}(1);
      t.month = date{1,1}(2);
      t.day = date{1,1}(3);
      t.hour=12; t.min=0; t.sec=0;
      t.UTC = obj.calib.model3D.UTC;
      sun = sun_position(t, obj.calib.model3D.Location);
      [x,y,z] = sph2cart(sun.azimuth*pi/180,(90-sun.zenith)*pi/180,1);
      p = [x y z]';
      pos = obj.camworld2im(R'*p);
end