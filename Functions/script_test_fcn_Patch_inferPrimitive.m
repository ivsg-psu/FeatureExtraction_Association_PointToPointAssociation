% script_test_fcn_Dataset_inferPrimitive.m
% This is a script to exercise the function:
% fcn_Dataset_inferPrimitive.m

% This script was written on 2022_01_21 by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history
%     2022_01_21
%     -- wrote the code

% Clean up the workspace
clearvars
clc

% Create the empty scalar structure for the patches
testPatches = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});

% Create a test patch that is orange and very vaguely round
testPatches(1).id = 'bafded';       % Give it a unique identifier
testPatches(1).color = [1 0.6 0.3]; % Set the color to a shade of orange
testPatches(1).primitive = 'irregular'; % By default, set the primitive shape to 'irregular'
% Add some points
testPatches(1).pointsX = [50.8640552995392; 58.2373271889401; 72.0622119815668; 68.6059907834101; 55.2419354838710];
testPatches(1).pointsY = [74.2335766423358; 66.9343065693431; 72.4817518248175; 85.3284671532847; 88.8321167883212];

% Create a test patch that is gray and very vaguely rectangular
testPatches(2,1).id = 'cefdab';       % Give it a unique identifier
testPatches(2).color = [0.6 0.6 0.6]; % Set the color to a shade of gray
testPatches(2).primitive = 'irregular'; % By default, set the primitive shape to 'irregular'
% Add some points
testPatches(2).pointsX = [69.7580645161290; 86.3479262672811; 83.1221198156682; 79.6658986175115; 72.5230414746544; 69.2972350230415];
testPatches(2).pointsY = [8.83211678832116; 20.2189781021898; 28.1021897810219; 24.8905109489051; 21.3868613138686; 17.8832116788321];

% Make some copies of the patch object to test different calls to
% inferPrimitive()
testPatches2 = testPatches;
testPatches3 = [testPatches; testPatches];

%% Find the best fit primitives to each patch, using different calls

% Find the best fit for all of the objects in the array of patches
fprintf('Before call A to inferPrimitive, testPatches.primitive and testPatches.primparams are:\n');
testPatches.primitive
testPatches.primparams
testPatches = fcn_Dataset_inferPrimitive(testPatches2);
fprintf('After call A to inferPrimitive, testPatches.primitive and testPatches.primparams are:\n');
testPatches.primitive
testPatches.primparams

% Create an empty figure ready for spatial plotting
figure(5)
clf         % Clear the figure for a new run of this script
hold on     % Turn on hold to plot patches, points, and fits
axis equal  % Set the axes to be equal spacing for spatial plotting
grid on     % Turn on the reference grid

% Plot the relevant patch before the fit, for reference
[h5,hpts5] = fcn_Dataset_plotPatch(testPatches2,5);    % Plot the patches with the definition points

% Find the best fit for all of the objects in the array of patches and plot
fprintf('Before call B to inferPrimitive, testPatches2.primitive and testPatches2.primparams are:\n');
testPatches2.primitive
testPatches2.primparams
testPatches2 = fcn_Dataset_inferPrimitive(testPatches2,5);
fprintf('After call B to inferPrimitive, testPatches2.primitive and testPatches2.primparams are:\n');
testPatches2.primitive
testPatches2.primparams

<<<<<<< HEAD:Functions/script_test_fcn_Patch_inferPrimitive.m
% Insert the points into the patch structure
testPatches(2) = fcn_Patch_insertPoints(testPatches(2),[x y]);

% Replot the data
delete(h(2))
delete(hpts(2))
[h(2),hpts(2)] = fcn_Patch_plotPatch(testPatches,5,2);
delete(hbb(2));
[testPatches,hbb(2)] = fcn_Patch_determineAABB(testPatches,5,2); % Determine and plot the axis-aligned bounding boxes
=======
% Find the best fit for the third object in the array of patches (without plotting)
fprintf('Before call C to inferPrimitive, testPatches3.primitive and testPatches3.primparams are:\n');
testPatches3.primitive
testPatches3.primparams
testPatches3(3) = fcn_Dataset_inferPrimitive(testPatches3(3));
fprintf('After call C to inferPrimitive, testPatches3.primitive and testPatches3.primparams are:\n');
testPatches3.primitive
testPatches3.primparams

% Create an empty figure ready for spatial plotting
figure(6)
clf         % Clear the figure for a new run of this script
hold on     % Turn on hold to plot patches, points, and fits
axis equal  % Set the axes to be equal spacing for spatial plotting
grid on     % Turn on the reference grid
>>>>>>> ce7f53c72d0990d0022a2d909815c1f656e5d143:Functions/script_test_fcn_Dataset_inferPrimitive.m

% Plot the relevant patch before the fit, for reference
[h6,hpts6] = fcn_Dataset_plotPatch(testPatches3(2),6);    % Plot the patches with the definition points

<<<<<<< HEAD:Functions/script_test_fcn_Patch_inferPrimitive.m
% Get rid of any existing primitives on the plot
if(exist('hprim','var'))
    delete(hprim)
end
% Fit and plot the best primitives (at the moment, this does only circles)
[testPatches, hprim] = fcn_Patch_inferPrimitive(testPatches,5);
=======
% Find the best fit for the second object in the array of patches and plot it
fprintf('Before call D to inferPrimitive, testPatches3.primitive and testPatches3.primparams are:\n');
testPatches3.primitive
testPatches3.primparams
testPatches3(2) = fcn_Dataset_inferPrimitive(testPatches3(2),6);
fprintf('After call D to inferPrimitive, testPatches3.primitive and testPatches3.primparams are:\n');
testPatches3.primitive
testPatches3.primparams
>>>>>>> ce7f53c72d0990d0022a2d909815c1f656e5d143:Functions/script_test_fcn_Dataset_inferPrimitive.m

