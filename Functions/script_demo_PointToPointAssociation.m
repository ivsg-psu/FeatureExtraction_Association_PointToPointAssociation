% This is a demonstration script to show the primary functionality of the
% point to point association library.

% Written
% 2023_05_23 

%% Prep the workspace
close all
clear all
clc

%% Show how input arguments are checked, fcn_Path_checkInputsToFunctions
% TO-DO - move debug tools into utilities and remove the checkinputs capability out of Path class

path_test = [4 1; 2 1];
fcn_Points_checkInputsToFunctions(path_test, 'path');

station_test = [1;3];
fcn_Points_checkInputsToFunctions(station_test, 'stations');

%% This is a script to exercise the function, fcn_Points_fillPointSampleSets

% Setup
% Call the function to fill in an array of "dataset" type
datasets_array = fcn_Points_fillPointSampleSets;

% Example 1
% We can save one of these as a single "dataset"
single_path = {datasets_array{1}};

fcn_Points_plotSetsXY(single_path,1);
text(single_path{1}(1,1),single_path{1}(1,2),'Start');

% Example 2
% illustrating the different XY points on different plots

for i_set = 1:length(datasets_array)
    fcn_Points_plotSetsXY(datasets_array(i_set),(i_set+1));
    text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
end

% Example 3
% All plots on the same graph

fcn_Points_plotSetsXY(datasets_array,5);

%% This is a script to exercise the function, fcn_Points_fillPointSetViaUserInputs

% Shows step by step instructions about using function

fig_num = 1;
h = figure(fig_num);
hold on;

num_iterations = input('How many point sets do you want to create? [Hit enter for default of 3]:','s');
if isempty(num_iterations)
    num_iterations = 3;
else
    num_iterations = str2double(num_iterations);
end
fprintf(1,'\n Filling in %.0d point sets.\n',num_iterations);
fprintf(1,'Instructions: \n');
fprintf(1,'Left click on the plot to create points. \n');
fprintf(1,'Right click on the plot to remove points \n');
fprintf(1,'Double click on the plot to end the set creation. \n');
fprintf(1,'When the last set is completed, another plot will be created to show results. \n');

    
% Initialize the paths_array
clear point_array
point_array{num_iterations} = [0 0];
for i_set = 1:num_iterations
    
    % Set the title header
    UserData.title_header = sprintf('Set %.0d of %.0d',i_set,num_iterations);
    
    % Save the results
    set(gcf,'UserData',UserData);
    
    pointsXY = fcn_Points_fillPointSetViaUserInputs(fig_num);
    point_array{i_set} = pointsXY;
end

clear data, close;
% Plot the results
fig_num = 13;
fcn_Points_plotSetsXY(point_array,fig_num);

%% This is a script to exercise the function, fcn_Points_plotSetsXY

% Example 1
% Plots XY positions of one dataset
fig_num = 1;

dataset = {};

dataset{1} = [1,21;
              2,2;
              3,3;
              4,1;
              5,4;
              6,32;
              7,0.5];

fcn_Points_plotSetsXY(dataset,fig_num)

% Example 2
% Plots XY positions of three datasets

dataset2 = {};
dataset2{1} = [1,21;
              2,2;
              3,3;
              4,1;
              5,4;
              6,32;
              7,0.5];

dataset2{2} = dataset{1}*1.8+0.3;

dataset2{3} = [9,2;
              8,0;
              2,5;
              1,1;
              6,4;
              10,31;
              2,1];

fcn_Points_plotSetsXY(dataset2,2)

% Example 3
% Plots from a sample data set inside fcn_Points_fillPointSampleSets

fig_num = 3;
datasets3 = fcn_Points_fillPointSampleSets;

fcn_Points_plotSetsXY(datasets3,fig_num);

for i_set = 1:length(datasets3)
    text(datasets3{i_set}(1,1),datasets3{i_set}(1,2),'Start');
end

%% This is a script to exercise the function, fcn_Points_pairXYdata

% Load up some data (simple xy points for now)
load testDatasetVehicle.mat

% Call the pairing function to obtain the matrix of paired XY data without
% a limiting radius
%[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyData{1},xyData{2});

% Define a maximum radius within which to consider data points as pairs
pairRadius = 1.0;

% Call the pairing function to obtain the matrix of paired XY data
[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyData{1},xyData{2},pairRadius);


% Now plot the matches with stars in various colors (to identify the pairs)
figure(17);
clf
hold on
quiver(pairedXYdata(:,1),pairedXYdata(:,2),pairedXYdata(:,3)-pairedXYdata(:,1),pairedXYdata(:,4)-pairedXYdata(:,2),0,'k','linewidth',1)

for i = 1:numMatches
    plot([pairedXYdata(i,1) pairedXYdata(i,3)],[pairedXYdata(i,2) pairedXYdata(i,4)],'*')
end
% Also plot the non-matches in data set A with circles in various colors
for i = 1+numMatches:nonMatchesA+numMatches
    plot(pairedXYdata(i,1),pairedXYdata(i,2),'o')
end
% Lastly, plot the non-matches in data set B with squares in various colors
for i = 1+numMatches+nonMatchesA:nonMatchesA+numMatches+nonMatchesA
    plot(pairedXYdata(i,3),pairedXYdata(i,4),'d')
end

%% This is a script to exercise the function, fcn_Points_calcPairStatistics

% Static test data for testing
pairedXYdata = [5.93317972350230,88.8321167883212,6.51881720430108,88.6921965317919;10.5414746543779,91.7518248175182,10.3718637992832,90.8598265895954;8.00691244239631,85.0364963503650,8.40053763440860,85.0794797687861;11.9239631336406,86.7883211678832,11.6263440860215,86.2355491329480;9.61981566820277,80.6569343065693,10.1030465949821,81.0332369942196;14.4585253456221,82.1167883211679,13.9560931899642,81.6112716763006;11.6935483870968,75.4014598540146,12.2535842293907,75.5419075144509;16.5322580645161,76.8613138686131,16.0170250896057,76.9869942196532;14.2281105990783,70.4379562043796,14.6729390681004,71.2066473988439;18.8364055299539,73.0656934306569,18.4363799283154,72.5072254335260;16.7626728110599,65.4744525547445,17.1818996415771,66.2933526011560;21.6013824884793,68.1021897810219,21.3037634408602,67.1604046242775;19.0668202764977,61.3868613138686,19.6908602150538,61.9580924855491;21.6013824884793,57.2992700729927,21.9310035842294,57.9118497109826;26.9009216589862,59.0510948905109,26.3216845878136,58.4898843930636;29.4354838709677,56.1313868613139,29.0098566308244,55.5997109826590;28.0529953917051,50.0000000000000,28.6514336917563,49.8193641618497;32.4308755760369,53.5036496350365,31.6980286738351,53.1430635838150;36.1175115207373,51.7518248175182,35.8198924731183,50.9754335260115;34.9654377880184,46.7883211678832,34.5654121863799,47.3627167630058;52.9377880184332,50.2919708029197,52.3073476702509,50.8309248554913;53.1682027649770,57.8832116788321,52.9345878136201,57.0447976878613;55.7027649769585,52.6277372262774,55.1747311827957,52.5650289017341;55.7027649769585,60.8029197080292,55.6227598566308,59.9349710982659;58.9285714285714,55.8394160583942,58.2213261648746,56.4667630057803;58.9285714285714,63.4306569343066,58.6693548387097,62.6806358381503;62.1543778801843,57.8832116788321,61.5367383512545,57.7673410404624;65.3801843317972,58.4671532846715,64.6729390681004,58.6343930635838;72.0622119815668,55.2554744525547,71.8413978494624,55.7442196531792;77.3617511520737,60.5109489051095,76.6801075268817,59.9349710982659;80.3571428571429,56.7153284671533,79.9059139784946,56.0332369942196;83.3525345622120,48.2481751824817,82.6836917562724,48.5187861271676;88.4216589861751,52.6277372262774,87.8808243727599,51.8424855491329;87.0391705069124,45.9124087591241,86.6263440860215,46.0621387283237;91.6474654377880,50.0000000000000,91.1066308243728,49.5303468208092;91.4170506912442,44.7445255474453,90.8378136200717,44.6170520231214;95.5645161290323,48.2481751824817,95.0492831541219,47.6517341040462;95.1036866359447,42.9927007299270,94.6012544802867,43.3164739884393];

% Run the statistics on the paired data
[errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(pairedXYdata);

% Post-process the shift to extract a distance and an angle for easy visual
% checking of the data
shiftDist = norm(meanShift,2);
shiftAngle = atan2(meanShift(2),meanShift(1));

%% This is a script to exercise the function, fcn_Points_addRadialNoise

% Example 1
% shows addition of radial noise

new_set1 = {[0.0 0.0]};
for i =  1:100
    datasets_out = fcn_Points_addRadialNoise(new_set1,0.5);
    fcn_Points_plotSetsXY(new_set1,1);
    fcn_Points_plotSetsXY(datasets_out,1);
    axis([-1 1 -1 1])
    axis square
end

% Example 2
% Addition of radial noise for one data set

new_set2 = {origXYdatasets{1}}; %#ok<CCAT1> 
datasets_out2 = fcn_Points_addRadialNoise(new_set2,1);
fcn_Points_plotSetsXY(new_set2,2);
fcn_Points_plotSetsXY(datasets_out2,2);


% Example 3
% Addition of radial noise for one data set

new_set3 = origXYdatasets;
datasets_out3 = fcn_Points_addRadialNoise(new_set3,1);
fcn_Points_plotSetsXY(new_set3,2);
fcn_Points_plotSetsXY(datasets_out3,2);

%% This is a script to exercise the function, fcn_Points_adjustPointSetStatistics 

% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;

%Corrupt one set of test data with only a systematic bias (for each)
biasedXYdataset = fcn_Points_adjustPointSetStatistics(origXYdatasets(1),[-0.1 0.1],zeros(1,2),zeros(1,2));

fh(1) = figure;
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')
    
horig(1) = fcn_Points_plotSetsXY(origXYdatasets(1),fh(1))
set(horig(1),'color','blue');
hcorr(1) = fcn_Points_plotSetsXY(biasedXYdataset,fh(1))
set(hcorr(1),'marker','*')
set(hcorr(1),'color','blue');
legend([horig(1) hcorr(1)],{'Original data','Biased data'})

% Call the pairing function to obtain pair the original data with the
% biased data
[pairedXYdataBias, numMatchesBias, nonMatchesABias, nonMatchesBBias] = fcn_Points_pairXYdata(origXYdatasets{1},biasedXYdataset{:});

% Calculate the statistics for the biased data set relative to the original
[errRMSBias,errVarBias,meanShiftBias] = fcn_Points_calcPairStatistics(pairedXYdataBias(1:numMatchesBias,:));

% Plot to provide a visual inspection of the bias
figure(17);
clf
hold on
quiver(pairedXYdataBias(:,1),pairedXYdataBias(:,2),pairedXYdataBias(:,3)-pairedXYdataBias(:,1),pairedXYdataBias(:,4)-pairedXYdataBias(:,2),0,'k','linewidth',1)
for i = 1:size(pairedXYdataBias,1)
    plot([pairedXYdataBias(i,1) pairedXYdataBias(i,3)],[pairedXYdataBias(i,2) pairedXYdataBias(i,4)],'*')
end

%% This is a script to exercise the function, fcn_Points_plotLaneMarkers

fig_num = 55;
path_test_laneMarkers = [4 1; 2 1; 1,21];


% path_test_laneMarkers = {[1,21;                                                 
%               2,2;
%               3,3;
%               4,1;
%               5,4;
%               6,32;
%               7,0.5]};                                                                   % -- Ask Dr. B

fcn_Points_plotLaneMarkers(path_test_laneMarkers,fig_num);

%% This is a script to exercise the function, fcn_Points_plotTrajectoryFromPath

fig_num = 555;
path_test_TrajfromPath = [4 1; 2 1; 1,21];


fcn_Points_plotTrajectoryFromPath(path_test_TrajfromPath, fig_num);                           % -- Ask Dr. B


%% This is a script to exercise the function, fcn_Patch_CalcCircularTrajectoryGeometry

% Vehicle trajectory information
vx = 20;        % longitudinal speed (m/s)
R = 20;         % path radius (m) with sign (+ left, - right)                                         -- Ask Dr. B
p0 = [10,-5];     % initial position of vehicle (m,m)
h0 = pi/2;     % initial heading of vehicle (rad)
a0 = 5*pi/180; % vehicle body slip angle (rad)
tf = 1;         % time horizon to check (s)
% Create a vector of the trajectory information
x0 = [p0'; h0; a0; vx; R];
% Vehicle dimensional information
vehicle.dr = 2.2;       % CG-front bumper distance (m)
vehicle.df = 3.5;        % CG-rear bumper distance (m)
vehicle.w = 2.3;        % vehicle width (m)

% Set up a new figure (or clear the existing one)
figure(1)
clf
hold on
grid on
axis equal

% Plot the vehicle CG at the initial location
plot(p0(1),p0(2),'ko','markersize',8)
% Plot a plus sign over the circle to make something akin to a CG symbol
plot(p0(1),p0(2),'k+','markersize',8)

% Calculate center point of vehicle trajectory circle
pc(1) = p0(1) + R*cos(h0+pi/2);
pc(2) = p0(2) + R*sin(h0+pi/2);

% Generate a trajectory for the given time horizon
del_t = 0.01;
t = (0:del_t:tf)';
theta = vx*t/R + h0-pi/2;
N = length(t);
pv = zeros(N,2);
pv(:,1) = R*cos(theta) + pc(1);
pv(:,2) = R*sin(theta) + pc(2);
% Plot the trajectory
plot(pv(:,1),pv(:,2),'k-.')

% Calculate the various pertinent radii and corner points of the vehicle
[radii,vehicleBB,radiiFlags] = fcn_Patch_CalcCircularTrajectoryGeometry(x0,vehicle);

% Parse out which radii are which
Rinside = sign(R)*radii(1);
RoutsideFront = sign(R)*radii(2);
RoutsideRear = sign(R)*radii(3);
Rmin = sign(R)*radii(6);
Rmax = sign(R)*radii(7);

% Calculate the offset to the angular position based on the back end of the
% vehicle (in order to plot the extra portion of the clearance curves)
rear_offset = 1.1*sign(R)*atan2(-vehicle.dr,sign(R)*R+vehicle.w/2);
% Create an angular range over which to plot the various radii
theta_arcs = linspace(rear_offset,theta(end),100);
% Plot the various radii 
plot(Rinside*cos(theta_arcs) + pc(1),Rinside*sin(theta_arcs) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(RoutsideFront*cos(theta_arcs) + pc(1),RoutsideFront*sin(theta_arcs) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(RoutsideRear*cos(theta_arcs) + pc(1),RoutsideRear*sin(theta_arcs) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(Rmin*cos(theta_arcs) + pc(1),Rmin*sin(theta_arcs) + pc(2),'r-');
plot(Rmax*cos(theta_arcs) + pc(1),Rmax*sin(theta_arcs) + pc(2),'b-');

% Finally, plot the vehicle body over the top of it all
plot(vehicleBB([1:4 1],1),vehicleBB([1:4 1],2),'k.-')
% If the inner tangent is the limiting inner radii, plot the tangent point
if 5 == radiiFlags(1)
    plot(vehicleBB(5,1),vehicleBB(5,2),'r*')
end

%% This is a script to exercise the function, fcn_Patch_checkCollisions

%-- Ask Dr. B

%% This is a script to exercise the function, fcn_Patch_checkStraightCollisions

%-- Ask Dr. B

%% This is a script to exercise the function, fcn_Patch_determineAABB

% Example 1
% Two test patches no overlap
% Create an empty figure ready for spatial plotting
figure(1)
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
testPatches(2,1).id = 'cefdab';
testPatches(2).color = [0.6 0.6 0.6];
testPatches(2).primitive = 'irregular';
testPatches(2).pointsX = [69.7580645161290; 86.3479262672811; 83.1221198156682; 79.6658986175115; 72.5230414746544; 69.2972350230415];
testPatches(2).pointsY = [8.83211678832116; 20.2189781021898; 28.1021897810219; 24.8905109489051; 21.3868613138686; 17.8832116788321];

testPatches = fcn_Patch_determineAABB(testPatches); % Determine the axis-aligned bounding boxes
[testPatches,hbb] = fcn_Patch_determineAABB(testPatches,1); % Determine and plot the axis-aligned bounding boxes

[h,hpts] = fcn_Patch_plotPatch(testPatches,1);    % Plot the patches with the definition points

%% This is a script to exercise the function, fcn_Patch_fillBarrelAtPoint

% Example 1
% Create and plot a barrel patch at [0,0] that is oriented at pi/8 angle

barrel1 = fcn_Patch_fillBarrelAtPoint([0 0],[0]);
grid on
fcn_Patch_plotPatch(barrel1,1);      % Plot orange patch with figure number 1
axis([-0.5 0.5 -0.5 0.5]);
axis square
% Example 2
% Create and plot a barrel patch at [0,0] that is oriented at pi/8 angle
barrel2 = fcn_Patch_fillBarrelAtPoint([1 1],[pi/8]);
fcn_Patch_plotPatch(barrel2,2);      % Plot orange patch with figure number 1
grid on
axis([0.5 1.5 0.5 1.5]);
axis square

% Example 3
% Create and plot 6 barrels

barrels = fcn_Patch_fillBarrelAtPoint([0 0; 1 1; 1.5 0; 0.5 0.7; 0 1.5;1.5 1.5],[pi/4; pi/16; pi/8; 0; pi/18; -pi]);
fcn_Patch_plotPatch(barrels,3);      % Plot orange patch with figure number 1
axis([-0.5 2 -0.5 2]);
axis square


%% This is a script to exercise the function, fcn_Patch_fillConeAtPoint

% Example 1
% Create and plot a cone patch at [0,0] that is oriented at pi/8 angle

cone1 = fcn_Patch_fillConeAtPoint([0 0],[0]);
grid on
fcn_Patch_plotPatch(cone1,1);      % Plot orange patch with figure number 1
axis([-2 2 -2 2]);
axis square
% Example 2
% Create and plot a cone patch at [0,0] that is oriented at pi/8 angle
cone2 = fcn_Patch_fillConeAtPoint([1 1],[pi/8]);
fcn_Patch_plotPatch(cone2,2);      % Plot orange patch with figure number 1
grid on
axis([-2 2 -2 2]);
axis square

% Example 3
% Create and plot 6 cones

cones = fcn_Patch_fillConeAtPoint([0 0; 1 1; 1.5 0; 0.5 0.7; 0 1.5;1.5 1.5],[pi/4; pi/16; pi/8; 0; pi/18; -pi]);
fcn_Patch_plotPatch(cones,3);      % Plot orange patch with figure number 1
axis([-0.5 2 -0.5 2]);
axis square

%% This is a script to exercise the function, fcn_Patch_fillPatchArrayViaUserInputs

% fcn_Patch_fillPatchArrayViaUserInputs
% A function for the user to click on the figure to generate XY data.
% Points are collected and plotted until the user double clicks. If the
% user right-clicks anywhere in the plot, the last point is deleted. Once
% the user double-clicks, the results are output from the function.

patchStruct = fcn_Patch_fillPatchArrayViaUserInputs(1);

%% This is a script to exercise the function, fcn_Patch_fillSamplePatches

% load twoPatches.mat
% samplePatches = fcn_Patch_fillSamplePatches(twoPatches.mat)          %        ---  Ask Dr. B

fcn_Patch_fillSamplePatches

%% This is a script to exercise the function, fcn_Patch_fillSquareAtPoint


square1 = fcn_Patch_fillSquareAtPoint([0 0],[0]);
grid on
fcn_Patch_plotPatch(square1,1);      % Plot orange patch with figure number 1
%axis([-2 2 -2 2]);
%axis square
% Example 2
% Create and plot a square patch at [0,0] that is oriented at pi/8 angle
square2 = fcn_Patch_fillSquareAtPoint([1 1],[pi/8]);
fcn_Patch_plotPatch(square2,2);      % Plot orange patch with figure number 1
grid on
%axis([-2 2 -2 2]);
%axis square

% Example 3
% Create and plot 6 squares

squares = fcn_Patch_fillSquareAtPoint([0 0; 1 1; 1.5 0; 0.5 0.7; 0 1.5;1.5 1.5],[pi/4; pi/16; pi/8; 0; pi/18; -pi]);
fcn_Patch_plotPatch(squares,3);      % Plot orange patch with figure number 1
%axis([-0.5 2 -0.5 2]);
%axis square

%% This is a script to exercise the function, fcn_Patch_fillType3BarricadeAtPoint

% Example 1
% Create and plot a Type3Barricade patch at [0,0] that is oriented at pi/8 angle

Type3Barricade1 = fcn_Patch_fillType3BarricadeAtPoint([0 0],[0]);
grid on
fcn_Patch_plotPatch(Type3Barricade1,1);      % Plot orange patch with figure number 1
axis equal

% Example 2
% Create and plot a Type3Barricade patch at [0,0] that is oriented at pi/8 angle
Type3Barricade2 = fcn_Patch_fillType3BarricadeAtPoint([1 1],[pi/8]);
fcn_Patch_plotPatch(Type3Barricade2,2);      % Plot orange patch with figure number 1
grid on
axis equal



% Example 3
% Create and plot 6 Type3Barricades

Type3Barricades = fcn_Patch_fillType3BarricadeAtPoint([0 0; 1 1; 1.5 0; 0.5 0.7; 0 1.5;1.5 1.5],[pi/4; pi/16; pi/8; 0; pi/18; -pi]);
fcn_Patch_plotPatch(Type3Barricades,3);      % Plot orange patch with figure number 1
axis auto
axis equal

%% This is a script to exercise the function, fcn_Patch_inferPrimitive

% -- Ask Dr. B

%% This is a script to exercise the function, fcn_Patch_insertPoints

% -- Ask Dr. B

%% This is a script to exercise the function, fcn_Patch_plotCircularCollisions



%% This is a script to exercise the function, fcn_Patch_plotCircularTrajectory

%% This is a script to exercise the function, fcn_Patch_plotPatch        
 
% Example 1
% Create and plot a test patch that is orange and very vaguely round


orangePatch = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});     % Ask Dr. B (What's h and What's hpts)

orangePatch(1).id = 'bafded';       % Give it a unique identifier
orangePatch(1).color = [1 0.6 0.3]; % Set the color to orange
orangePatch(1).primitive = 'irregular'; % By default, set the primitive shape to 'irregular'
% Add some points
orangePatch(1).pointsX = [50.8640552995392; 58.2373271889401; 72.0622119815668; 68.6059907834101; 55.2419354838710];
orangePatch(1).pointsY = [74.2335766423358; 66.9343065693431; 72.4817518248175; 85.3284671532847; 88.8321167883212];

[h1, hpts1] = fcn_Patch_plotPatch(orangePatch,1);      % Plot orange patch with figure number 1
f1 = gcf;
fprintf("\nCall A: At this point, you should have one figure open with an orange patch plotted and\nthe 1x1 handle variables h1 and hpts1 in the workspace.\n\n");

% Example 2
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


% Example 3
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

% Plot all of the patches (with points) in the patch array and save a vector
% of handles to each of the patch objects and point plot objects
[h3, hpts3] = fcn_Patch_plotPatch(testPatches,4);
f4 = gcf;
fprintf("\nCall D: At this point, you should have one figure open with an orange, a gray, and a blue patch plotted\nand the 3x1 handle vectors h2 and hpts2 in the workspace.\n\n");
whos

%% This is a script to exercise the function, fcn_Patch_plotStraightCollisions

%% This is a script to exercise the function, fcn_Patch_plotStraightTrajectory






