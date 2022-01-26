% script_test_fcn_Dataset_plotPatch.m
% This is a script to exercise the function:
% fcn_Dataset_plotPatch.m

% This script was written on 2022_01_26 by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history
%     2022_01_26
%     -- wrote the code

% Clean up the workspace
close all
clearvars
clc

pause(0.2);

% Create the empty scalar structure for the patches
testPatches = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});

% Create a test patch that is orange and very vaguely round
testPatches(1).id = 'bafded';       % Give it a unique identifier
testPatches(1).color = [1 0.6 0.3]; % Set the color to orange
testPatches(1).primitive = 'irregular'; % By default, set the primitive shape to 'irregular'
% Add some points
testPatches(1).pointsX = [50.8640552995392; 58.2373271889401; 72.0622119815668; 68.6059907834101; 55.2419354838710];
testPatches(1).pointsY = [74.2335766423358; 66.9343065693431; 72.4817518248175; 85.3284671532847; 88.8321167883212];

% Create a test patch that is gray and very vaguely rectangular
testPatches(2,1).id = 'cefdab';
testPatches(2).color = [0.6 0.6 0.6];
testPatches(2).primitive = 'irregular';
testPatches(2).pointsX = [69.7580645161290; 86.3479262672811; 83.1221198156682; 79.6658986175115; 72.5230414746544; 69.2972350230415];
testPatches(2).pointsY = [8.83211678832116; 20.2189781021898; 28.1021897810219; 24.8905109489051; 21.3868613138686; 17.8832116788321];

%% Perform the patch plotting

% Plot all of the patches (with points) in the patch array on a new figure
% with no return arguments
fcn_Dataset_plotPatch(testPatches);
f1 = gcf;
fprintf("\nCall A: At this point, you should have one figure open with an orange and a gray patch plotted.\n\n");
pause();

% Plot all of the patches (with points) in the patch array and save a vector
% of handles to each of the patch objects and point plot objects
[h2,hpts2] = fcn_Dataset_plotPatch(testPatches);
f2 = gcf;
fprintf("\nCall B: At this point, you should have a second figure open with an orange and a gray patch plotted\nand the 2x1 handle vectors h2 and hpts2 in the workspace.\n\n");
whos
pause();

% Plot a single patch (with points) from the patch array and save the
% vector of handles to the patch object and point plot object
[h3,hpts3] = fcn_Dataset_plotPatch(testPatches(2));
f3 = gcf;
fprintf("\nCall C: At this point, you should have a third figure open with only the gray patch plotted and\nthe 1x1 handle variables h3 and hpts3 in the workspace.\n\n");
whos
pause();

%% Create an empty figure ready for spatial plotting
figure(4)
clf         % Clear the figure for a new run of this script
axis equal  % Set the axes to be equal spacing for spatial plotting
grid on     % Turn on the reference grid

% Plot all of the patches (with points) in the patch array into the
% provided figure and save the vectors of handles
[h4,hpts4] = fcn_Dataset_plotPatch(testPatches,4);
f4 = gcf;
fprintf("\nCall D: At this point, you should have figure 4 (specifically) open with an orange and a gray patch plotted.\n\n");
pause();

% Add the first patch on figure 3 (which previously only had the plot of
% the second patch) and add the handles for the plot objects to the handle
% vectors
figure(f3)          % Foreground the third figure for the user to see
[h3(2,1),hpts3(2,1)] = fcn_Dataset_plotPatch(testPatches,f3,1);
fprintf("\nCall E: At this point, the third figure created (initially with only the gray patch) should have the orange patch added.\nThe handle variables h3 and hpts3 should now be 2x1 to reflect the new plot objects.\n\n");
whos
pause();

% Delete the first object from the second figure created
figure(f2)          % Foreground the second figure for the user to see
delete(h2(1))       % Delete the patch object
delete(hpts2(1))    % Delete the points on the object
h2 = h2(2);         % Remove the patch handle from the patch handle vector
hpts2 = hpts2(2);   % Remove the points handle from the points handle vector
fprintf("\nCall F: This is not a call to the function, per se, but on the second figure, the orange patch should be gone.\nThe object handles h2 and hpts2 should now be 1x1 to reflect this.\n\n");
whos
pause();