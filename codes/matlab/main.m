% close all; 
clear all; 
clc;

%% settings
prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj\';
% prj_path = 'E:\ms_proj\';
addpath(genpath([prj_path,'\old_codes\vismolib_v6']));
addpath(genpath([prj_path,'\lib']));
addpath(prj_path);


date = '2015_08_03';
s = vml('cavriglia',date);
