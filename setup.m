%% ^_^ Welcome to Tensegrity Finite Element Method(TsgFEM) software! ^_^ %%
% SETUP file to be run only the first time
%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

%% Add all necessary functions to MATLAB path%% set up the workspace
clear all;
close all;
clc;
% add the function libraries
addpath( genpath('codes') );
% add the Software Verification and Examples
addpath( genpath('examples') );

addpath( genpath('function') );
