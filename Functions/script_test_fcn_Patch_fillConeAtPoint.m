% script_test_fcn_Patch_fillConeAtPoint.m
% This is a script to exercise the function:
% fcn_Patch_fillConeAtPoint.m

% This script was written on 2022_07_18 by Shashank Garikipati


%% Revision history
% 2022_07_18 - Shashank Garikipati
% -- wrote the code


% Clean up the workspace

clc, close all

% script_test_fcn_Patch_fillConeAtPoint.m
% This is a script to exercise the function:
% fcn_Patch_fillConeAtPoint.m

% This script was written on 2022_07_05 by Shashank Garikipati



%% Example 1
% Create and plot a cone patch at [0,0] that is oriented at pi/8 angle

cone1 = fcn_Patch_fillConeAtPoint([0 0],[0]);
grid on
fcn_Patch_plotPatch(cone1,1);      % Plot orange patch with figure number 1
axis([-2 2 -2 2]);
axis square
%% Example 2
% Create and plot a cone patch at [0,0] that is oriented at pi/8 angle
cone2 = fcn_Patch_fillConeAtPoint([1 1],[pi/8]);
fcn_Patch_plotPatch(cone2,2);      % Plot orange patch with figure number 1
grid on
axis([-2 2 -2 2]);
axis square

%% Example 3
% Create and plot 6 cones

cones = fcn_Patch_fillConeAtPoint([0 0; 1 1; 1.5 0; 0.5 0.7; 0 1.5;1.5 1.5],[pi/4; pi/16; pi/8; 0; pi/18; -pi]);
fcn_Patch_plotPatch(cones,3);      % Plot orange patch with figure number 1
axis([-0.5 2 -0.5 2]);
axis square