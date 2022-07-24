% script_test_fcn_Patch_fillType3BarricadeAtPoint.m
% This is a script to exercise the function:
% fcn_Patch_fillType3BarricadeAtPoint.m

% This script was written on 2022_07_05 by Shashank Garikipati


%% Revision history
% 2022_01_26 - Shashank Garikipati
% -- wrote the code
% 2022_07_18 - Shashank Garikipati
% Cahnged functionality of code from gi input to matrix input

% Clean up the workspace

clc, close all

%% Example 1
% Create and plot a Type3Barricade patch at [0,0] that is oriented at pi/8 angle

Type3Barricade1 = fcn_Patch_fillType3BarricadeAtPoint([0 0],[0]);
grid on
fcn_Patch_plotPatch(Type3Barricade1,1);      % Plot orange patch with figure number 1
axis equal

%% Example 2
% Create and plot a Type3Barricade patch at [0,0] that is oriented at pi/8 angle
Type3Barricade2 = fcn_Patch_fillType3BarricadeAtPoint([1 1],[pi/8]);
fcn_Patch_plotPatch(Type3Barricade2,2);      % Plot orange patch with figure number 1
grid on
axis equal



%% Example 3
% Create and plot 6 Type3Barricades

Type3Barricades = fcn_Patch_fillType3BarricadeAtPoint([0 0; 1 1; 1.5 0; 0.5 0.7; 0 1.5;1.5 1.5],[pi/4; pi/16; pi/8; 0; pi/18; -pi]);
fcn_Patch_plotPatch(Type3Barricades,3);      % Plot orange patch with figure number 1
axis auto
axis equal






