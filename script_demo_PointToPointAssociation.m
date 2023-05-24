% This is a demonstration script to show the primary functionality of the
% point to point association library.

% Written
% 2023_05_23 A. Batchu, aneeshb@psu.edu

%% Prep the workspace
close all
clc

%% Show how input arguments are checked, fcn_Path_checkInputsToFunctions
% TO-DO - move debug tools into utilities and remove the checkinputs capability out of Path class

path_test = [4 1; 2 1];
fcn_Points_checkInputsToFunctions(path_test, 'path');