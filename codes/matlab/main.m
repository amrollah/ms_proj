% close all; 
clear all; 
clc;

%% settings
% prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj';
prj_path = 'E:\ms_proj\';
addpath(genpath([prj_path,'\old_codes\vismolib_v5']));
addpath(genpath([prj_path,'\lib']));
addpath(prj_path);


date = '201_02_01';
s = vmlSeq('cavriglia',date);
