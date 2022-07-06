% script_test_fcn_Patch_fillBarrelAtPoint.m
% This is a script to exercise the function:
% fcn_Patch_fillBarrelAtPoint.m

% This script was written on 2022_07_05 by Shashank Garikipati


%% Revision history
% 2022_01_26 - Shashank Garikipati
% -- wrote the code

% Clean up the workspace

clc, close all

%% Example 1
% Create and plot a barrel patch at [2,3] that is oriented at pi/8 angle
figure(1);
axis([-5 5 -5 5]);
barrel1 = fcn_Patch_fillBarrelAtPoint(2,pi/8);
barrel2 = fcn_Patch_fillBarrelAtPoint(2,0);
barrel3 = fcn_Patch_fillBarrelAtPoint(2,pi/16)

fcn_Patch_plotPatch(barrel1,1);      % Plot orange patch with figure number 1
fcn_Patch_plotPatch(barrel2,1);      % Plot orange patch with figure number 1
fcn_Patch_plotPatch(barrel3,1);      % Plot orange patch with figure number 1

