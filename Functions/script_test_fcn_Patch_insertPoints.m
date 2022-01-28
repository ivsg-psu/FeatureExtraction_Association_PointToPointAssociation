% script_test_fcn_Patch_insertPoints.m
% This is a script to exercise the function:
% fcn_Patch_insertPoints.m

% This script was written on 2022_01_26 by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history
%     2022_01_26
%     -- wrote the code

% Clean up the workspace
clearvars

% Create an empty figure ready for spatial plotting
figure(5)
clf         % Clear the figure for a new run of this script
axis equal  % Set the axes to be equal spacing for spatial plotting
grid on     % Turn on the reference grid

% Load some sample patch data to work on
testPatches = fcn_Patch_fillSamplePatches;

% Plot the patches with the vertices for visualization
[h,hpts] = fcn_Patch_plotPatch(testPatches,5);    


%% Insert new points into the patch object

% Create test points using the current figure. Click to add points, hit
% enter to finish.
[x,y] = ginput;

% Insert the points into the patch structure. Save the returned patch
% structure to the same location to overwrite the old data.
testPatches(1) = fcn_Patch_insertPoints(testPatches(1),[x y]);

% Replot the data by deleting the graphics objects referenced by the first
% elements in the handle vectors, matching with the first patch object in
% the structure array. (If the plotting order was different, the delete
% and reassignment of the handles would need to match the plotting order.)
delete(h(1))
delete(hpts(1))
[h(1),hpts(1)] = fcn_Patch_plotPatch(testPatches,5,1);