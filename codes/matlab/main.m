close all; 
clear all; 
clc;

%% settings
prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj';
prj_path = 'E:\ms_proj';
addpath(genpath([prj_path,'\old_codes\vismolib_v4']));
addpath(genpath([prj_path,'\lib']));
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');

%
date = '2015_08_03';
s = vmlSeq(date);