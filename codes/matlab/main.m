% close all; 
clear all; 
clc;

%% settings
% prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj';
prj_path = 'E:\ms_proj\';
addpath(genpath([prj_path,'\old_codes\vismolib_v5']));
addpath(genpath([prj_path,'\lib']));
% run([prj_path 'local_conf.m']);
% conf = evalin('base','VMLCONF');


%
date = '2015_12_19';
s = vmlSeq('cavriglia',date);
