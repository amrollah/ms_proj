function pscp(host,files,pw)
if nargin<2, error('Not enough input arguments'); end;
if nargin<3, pw = input('Password: ','s'); end;
clc;
if isempty(findstr(files,'/')), files = [files '/m*.mat']; end;
str = ['pscp -pw ' pw ' poland@' host ':/usr/local/users/poland/matlab/' files ' .'];
str1 = ['pscp -pw *** poland@' host ':/usr/local/users/poland/matlab/' files ' .'];
disp(str1);
system(str);
