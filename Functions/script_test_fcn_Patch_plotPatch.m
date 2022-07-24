% script_test_fcn_Patch_plotPatch.m
% This is a script to exercise the function:
% fcn_Patch_plotPatch.m

% This script was written on 2022_01_26 by C. Beal
% Questions or comments? cbeal@bucknell.edu


%% Revision history
% 2022_01_26 - C. Beal
% -- wrote the code
% 2022_07_03 - Shashank Garikipati
% -- organized script into 4 examples and added patch for demonstration

% Clean up the workspace

clc
clear vars
close

%% Example 1
% Create and plot a test patch that is orange and very vaguely round

orangePatch = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});

orangePatch(1).id = 'bafded';       % Give it a unique identifier
orangePatch(1).color = [1 0.6 0.3]; % Set the color to orange
orangePatch(1).primitive = 'irregular'; % By default, set the primitive shape to 'irregular'
% Add some points
orangePatch(1).pointsX = [50.8640552995392; 58.2373271889401; 72.0622119815668; 68.6059907834101; 55.2419354838710];
orangePatch(1).pointsY = [74.2335766423358; 66.9343065693431; 72.4817518248175; 85.3284671532847; 88.8321167883212];

[h1, hpts1] = fcn_Patch_plotPatch(orangePatch,1);      % Plot orange patch with figure number 1
f1 = gcf;
fprintf("\nCall A: At this point, you should have one figure open with an orange patch plotted and\nthe 1x1 handle variables h1 and hpts1 in the workspace.\n\n");
pause();

%% Example 2
% Create and plot a test patch that is gray and very vaguely rectangular

fig_num = 2;

grayPatch = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});

% Create a test patch that is gray and very vaguely rectangular
grayPatch(1).id = 'cefdab';
grayPatch(1).color = [0.6 0.6 0.6];
grayPatch(1).primitive = 'irregular';
grayPatch(1).pointsX = [69.7580645161290; 86.3479262672811; 83.1221198156682; 79.6658986175115; 72.5230414746544; 69.2972350230415];
grayPatch(1).pointsY = [8.83211678832116; 20.2189781021898; 28.1021897810219; 24.8905109489051; 21.3868613138686; 17.8832116788321];


[h2, hpts2] = fcn_Patch_plotPatch(grayPatch,fig_num);      % Plot gray patch with figure number
f2 = gcf;
fprintf("\nCall B: At this point, you should have one figure open with an gray patch plotted and\nthe 1x1 handle variables h2 and hpts2 in the workspace.\n\n");
pause();

%% Example 3
% Save orange, gray, and new blue triangle patch into one structure and plot them on the same plot

testPatches = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});

testPatches(1,:) = orangePatch;    % Run examples 1 and 2 before runnning example 3
testPatches(2,:) = grayPatch;

% Create a new blue triangle patch into the third row of the testPatches
% struct array
testPatches(3,1).id = 'dwqbed';       % Give it a unique identifier
testPatches(3).color = [0.2 0.3 0.8]; % Set the color to blue
testPatches(3).primitive = 'irregular'; % By default, set the primitive shape to 'irregular'
% Add some points
testPatches(3).pointsX = [55; 60; 65];
testPatches(3).pointsY = [30; 50; 20];


% Plot all of the patches (with points) in the patch array on a new figure
% with no return arguments
fcn_Patch_plotPatch(testPatches,3);
f3 = gcf;
fprintf("\nCall C: At this point, you should have one figure open with an orange, a gray, and a blue patch plotted.");
whos
pause();

% Plot all of the patches (with points) in the patch array and save a vector
% of handles to each of the patch objects and point plot objects
[h3, hpts3] = fcn_Patch_plotPatch(testPatches,4);
f4 = gcf;
fprintf("\nCall D: At this point, you should have one figure open with an orange, a gray, and a blue patch plotted\nand the 3x1 handle vectors h2 and hpts2 in the workspace.\n\n");
whos
pause();


%% Example 4
% Create an empty figure ready for spatial plotting
figure(5)
clf         % Clear the figure for a new run of this script
axis equal  % Set the axes to be equal spacing for spatial plotting
grid on     % Turn on the reference grid

% Plot all of the patches (with points) in the patch array into the
% provided figure and save the vectors of handles
[h4,hpts4] = fcn_Patch_plotPatch(testPatches,4);
f5 = gcf;
fprintf("\nCall E: At this point, you should have figure 4 (specifically) open with an orange, gray, and blue patch plotted.\n\n");
pause();

%% Example 5
% Add the first patch on figure 2 (which previously only had the plot of
% the second patch) and add the handles for the plot objects to the handle
% vectors
figure(f2)          % Foreground the second figure for the user to see
[h2(2,1),hpts2(2,1)] = fcn_Patch_plotPatch(testPatches,f2,1);
fprintf("\nCall F: At this point, the third figure created (initially with only the gray patch) should have the orange patch added.\nThe handle variables h2 and hpts2 should now be 2x1 to reflect the new plot objects.\n\n");
whos
pause();


%% Example 6
% Delete the first object from the second figure created
figure(f2)          % Foreground the third figure for the user to see
delete(h2(1))       % Delete the patch object
delete(hpts2(1))    % Delete the points on the object
h3 = h3(2);         % Remove the patch handle from the patch handle vector
hpts3 = hpts2(2);   % Remove the points handle from the points handle vector
fprintf("\nCall G: This is not a call to the function, per se, but on the second figure, the orange patch should be gone.\nThe object handles h2 and hpts2 should now be 1x1 to reflect this.\n\n");
whos
pause();