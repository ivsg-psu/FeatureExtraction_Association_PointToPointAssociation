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
testPatches(2).id = 'cefdab';
testPatches(2).color = [0.6 0.6 0.6];
testPatches(2).primitive = 'irregular';
testPatches(2).pointsX = [69.7580645161290; 86.3479262672811; 83.1221198156682; 79.6658986175115; 72.5230414746544; 69.2972350230415];
testPatches(2).pointsY = [8.83211678832116; 20.2189781021898; 28.1021897810219; 24.8905109489051; 21.3868613138686; 17.8832116788321];

%testPatches = fcn_Patch_determineAABB(testPatches); % Determine the axis-aligned bounding boxes
[testPatches,hbb] = fcn_Patch_determineAABB(testPatches,5); % Determine and plot the axis-aligned bounding boxes

[h,hpts] = fcn_Patch_plotPatch(testPatches,5);    % Plot the patches with the definition points


%% Insert a point into the patch object

% Create test points
[x,y] = ginput;

% Insert the points into the patch structure
testPatches(2) = fcn_Patch_insertPoints(testPatches(2),[x y]);

% Replot the data
delete(h(2))
delete(hpts(2))
[h(2),hpts(2)] = fcn_Patch_plotPatch(testPatches,5,2);
delete(hbb(2));
%[testPatches,hbb(2)] = fcn_Patch_determineAABB(testPatches,5,2); % Determine and plot the axis-aligned bounding boxes
