%% Documentation   
% Standalone script to post-process linearization files from OpenFAST for one operating point.
%
%% Initialization
clear all; close all; clc;
restoredefaultpath;
addpath(genpath('C:/Work/FAST/matlab-toolbox/')); % TODO adapt me

%% Script Parameters
basename     = './Main'; % Basename, lin files will be assumed to be base_name.i.lin
BladeLen     = 61.5    ; % Blade length, used to tune relative modal energy [m]
TowerLen     = 87.6    ; % Tower length, used to tune relative modal energy [m]
nModesMax    = 10      ; % Maximum number of modes to be shown
nCharMaxDesc = 50      ; % Maximum number of characters for description written to screen

modeVizFileName='./Main.ModeShapeVTK.postMBC';

%% Derived parameters
lin_files = cell(3,1);
ID = [1,12,24];
for i = 1:length(ID)
    lin_files{i} = sprintf('%s.%d.lin', basename, ID(i));
end
lin_files

%% Performing MBC (NOTE: not stricly necessary without rotation)
[mbc_data, matData, FAST_linData] = fx_mbc3( lin_files , modeVizFileName);
[CampbellData] = campbell_diagram_data(mbc_data, BladeLen, TowerLen); %, xlsFileName);

%% Outputs to screen and save to file
[Freq, Damp] = printCampbellDataOP(CampbellData);

