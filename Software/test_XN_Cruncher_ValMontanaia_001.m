% Made by D. Zuliani 2025/04/12
% This is a script to test XN_Cruncher on the dataset "Val Montanaia 001"
% included in the "Data" dir of this distribution. This dataset format is
% ASCII.
%
%% CLEANING
clear all;
close all;
format long g;
%
%% Setting SLASH for computer dependent PATHS
if ispc
    SLASH_TYPE = '\';
else
    SLASH_TYPE = '/';
end
%
%% Setting script path
[SCRIPTPATH, ~, ~] = fileparts(mfilename('fullpath'));
DATAPATH_IN     = [SCRIPTPATH,SLASH_TYPE,'..',SLASH_TYPE,'DataIn'];
DATAPATH_OUT    = [SCRIPTPATH,SLASH_TYPE,'..',SLASH_TYPE,'DataOut'];
%
%% FILELIST
% FILELIST is a structure array which includes 3 full filename with pathname.
% Each file must contain a component of a velocimeter sensor. The list must be
% sorted as below:
% 1st component: East - West
% 2nd component: North - South
% 3th component: Vertical
% The data format allowed are both ascii and sac. 
% e.g. 
% FILELIST = {'/home/myuser/EHE.vel','/home/myuser/EHN.vel','/home/myuser/EHZ.vel'};
% FILELIST = {'/Users/myuser/EHE.vel','/Users/myuser/EHN.vel','/Users/myuser/EHZ.vel'};
% FILELIST = {'C:\Users\dzuliani\EHE.vel','C:\Users\dzuliani\EHN.vel','C:\Users\dzuliani\EHZ.vel'};
FILELIST = {[DATAPATH_IN,SLASH_TYPE,'ValMontanaia_001_E-W.trc'],...
    [DATAPATH_IN,SLASH_TYPE,'ValMontanaia_001_N-S.trc'],...
    [DATAPATH_IN,SLASH_TYPE,'ValMontanaia_001_Z.trc']};
FILE_MATLAB_OUT = [DATAPATH_OUT,SLASH_TYPE,'ValMontanaia_001.mat'];
%
%% Running XN_Cruncher
XN_DATA=XN_Cruncher(FILELIST,FILE_MATLAB_OUT);