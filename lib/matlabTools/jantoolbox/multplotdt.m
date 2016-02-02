function [haxes, hlines] = multplotdt(tt,yy,linestyles,titles,varargin)
if nargin<3, linestyles=''; end
if nargin<4, titles={}; end
[haxes, hlines] = multplot(tt,yy,linestyles,titles,varargin{:});
dynamicDateTicks(haxes,'link');
