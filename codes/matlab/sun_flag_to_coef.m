function coefs = sun_flag_to_coef(sun_flag)
  if ~isrow(sun_flag)
      sun_flag = (sun_flag);
  end
   coefs = arrayfun(@(flag) convert(flag),sun_flag);
    function coef = convert(sun_flag)
      switch lower(sun_flag)
        case 1 % no visible sun, DNI~0
          coef = 0;
        case 4 % complete star shape sun, DNI~1
          coef = 1;
        otherwise
          coef = NaN;
      end
    end
 end