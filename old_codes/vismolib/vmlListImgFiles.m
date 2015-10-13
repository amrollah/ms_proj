function [files,folder,tt] = vmlListImgFiles(folder,hours,conf)
% [files,dt] = vmlListImgFiles(folder,[hours],[basefolder])
% List camera images in specified folder.
% If hours is specified, then only file names between hours(1) and
% hours(2) will be returned, including some lead minutes.
% If basefolder is unspecified, it will be taken from VMLCONF.
% The function returns the file names, folder, and the numeric date-times.
%
% 2015-04-07, Jan Poland, ABB.CH-RD.C1

if nargin<2, hours = []; end
if nargin<3, conf = evalin('base','VMLCONF'); end
if isempty(hours), hours = [0 24]; end

folder = [conf.basefolder folder '\'];
files = dir([folder '*' conf.img_ext]);
if isempty(files), error(['no data in folder ' folder]); end
files = cellfun(@(x)strrep(x,conf.img_ext,''),{files.name},'UniformOutput',false');
t = cellfun(@(s)textscan(s,'%s','delimiter','_'),files);
N = sum(~isnan(str2double(t{1})));
t=cellfun(@(x)cellfun(@str2double,x(1:N)),t,'UniformOutput',0);
t = [t{:}];
j = t(4,:)<hours(2) & (t(4,:)>=hours(1) | (t(4,:)==hours(1)-1 & t(5,:)>=60-conf.pad_minutes));
t = t(:,j);
files = files(j);
tt=datenum(t(1:6,:)')';
if N>=7, tt = tt+t(7,:)/86400000; end
