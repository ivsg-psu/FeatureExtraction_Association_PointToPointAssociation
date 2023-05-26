%% Introduction to and Purpose of the PointToPoint Code
% This is a demonstration script to show the primary functionality of the
% point to point association library.
%
% This is the explanation of the code that can be found by running
%       script_demo_PointToPointAssociation.m
%
% This is a script to demonstrate the functions within the PointToPoint code
% library. This code repo is typically located at:
%   https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation
%
% If you have questions or comments, please contact Sean Brennan at
% sbrennan@psu.edu or Aneesh Batchu, aneeshb@psu.edu

% TO-DO - move debug tools into utilities and remove the checkinputs capability out of Path class


%% Revision History:
% 2023_05_23 A. Batchu, aneeshb@psu.edu
% -- Started the master script
% 2023_05_24: S. Brennan, sbrennan@psu.edu
% -- added automatic loading of utilities

%% Prep the workspace
close all
clc

%% Dependencies and Setup of the Code
% The code requires several other libraries to work, namely the following
% 
% * DebugTools - the repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools
% * PathClassLibrary - the repo can be found at: https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary
% * GPS - this is the library that converts from ENU to/from LLA
% * GetUserInputs - this 
% List what libraries we need, and where to find the codes for each
clear library_name library_folders library_url

ith_library = 1;
library_name{ith_library}    = 'DebugTools_v2023_04_22';
library_folders{ith_library} = {'Functions','Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/archive/refs/tags/DebugTools_v2023_04_22.zip';

ith_library = ith_library+1;
library_name{ith_library}    = 'PathClass_v2023_02_01';
library_folders{ith_library} = {'Functions'};                                
library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/blob/main/Releases/PathClass_v2023_02_01.zip?raw=true';

%% Clear paths and folders, if needed
if 1==0

   fcn_INTERNAL_clearUtilitiesFromPathAndFolders;

end

%% Do we need to set up the work space?
if ~exist('flag_PointToPoint_Folders_Initialized','var')
    this_project_folders = {'Functions'}; % {'Functions','Data'};
    fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders);  
    flag_PointToPoint_Folders_Initialized = 1;
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Getting%20Started
% 
%    _____      _   _   _                _____ _             _           _ 
%   / ____|    | | | | (_)              / ____| |           | |         | |
%  | |  __  ___| |_| |_ _ _ __   __ _  | (___ | |_ __ _ _ __| |_ ___  __| |
%  | | |_ |/ _ \ __| __| | '_ \ / _` |  \___ \| __/ _` | '__| __/ _ \/ _` |
%  | |__| |  __/ |_| |_| | | | | (_| |  ____) | || (_| | |  | ||  __/ (_| |
%   \_____|\___|\__|\__|_|_| |_|\__, | |_____/ \__\__,_|_|   \__\___|\__,_|
%                                __/ |                                     
%                               |___/                                      
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This is a script to exercise the function, fcn_Points_plotSetsXY


% Example 1
% Plots XY positions of one dataset
fig_num = 8;

dataset = {};

dataset{1} = [1,21;
              2,2;
              3,3;
              4,1;
              5,4;
              6,32;
              7,0.5];

fcn_Points_plotSetsXY(dataset,fig_num)
title('fcn_Points_plotSetsXY: Example plotting XY positions of one dataset','Interpreter','none')

% Example 2
% Plots XY positions of three datasets

fig_num = 9;

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

fcn_Points_plotSetsXY(dataset2,fig_num)
title('fcn_Points_plotSetsXY: Example plotting XY positions of three datasets','Interpreter','none')

% Example 3
% Plots from a sample data set inside fcn_Points_fillPointSampleSets

fig_num = 10;
datasets3 = fcn_Points_fillPointSampleSets;

fcn_Points_plotSetsXY(datasets3,fig_num);
title('fcn_Points_plotSetsXY: Example plotting sample data set inside fcn_Points_fillPointSampleSets','Interpreter','none')
for i_set = 1:length(datasets3)
    text(datasets3{i_set}(1,1),datasets3{i_set}(1,2),'Start');
end

%% This is a script to exercise the function, fcn_Points_fillPointSampleSets

% Setup
% Call the function to fill in an array of "dataset" type
datasets_array = fcn_Points_fillPointSampleSets;


% Example 1
% We can save one of these as a single "dataset"
fig_num = 1;
single_path = {datasets_array{1}};

fcn_Points_plotSetsXY(single_path,fig_num);
text(single_path{1}(1,1),single_path{1}(1,2),'Start');
title('fcn_Points_plotSetsXY: Example plotting of one dataset ','Interpreter','none')  
% Example 2
% illustrating the different XY points on different plots

for i_set = 1:length(datasets_array)
    fig_num = i_set+1;
    fcn_Points_plotSetsXY(datasets_array(i_set),fig_num);
    text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
    title(sprintf('fcn_Points_plotSetsXY: Example plotting of dataset: %.0d',i_set),'Interpreter','none'); 
end

% Example 3
% All plots on the same graph
fig_num = 5;
fcn_Points_plotSetsXY(datasets_array,fig_num);
title('fcn_Points_plotSetsXY: Example plotting of all datasets in one call','Interpreter','none')



%% This is a script to exercise the function, fcn_Points_fillPointSetViaUserInputs

% Shows step by step instructions about using function

fig_num = 6;
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
    title('fcn_Points_fillPointSetViaUserInputs: Plotting the points via user inputs','Interpreter','none')
    point_array{i_set} = pointsXY;
end

clear data, close;
% Plot the results
fig_num = 7;
fcn_Points_plotSetsXY(point_array,fig_num);
title('Plotting the point array (points plotted via user inputs) in a new figure','Interpreter','none')


%% This is a script to exercise the function, fcn_Points_adjustPointSetStatistics 

% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;
biasx = -2.2;
biasy = 1.2;

noiseMean = zeros(1,2);
noiseVariance = zeros(1,2);
%Corrupt one set of test data with only a systematic bias (for each)
biasedXYdataset = fcn_Points_adjustPointSetStatistics(origXYdatasets(1),[biasx biasy],noiseMean,noiseVariance);

% Plotting the original data and biased data
% fh(1) = figure;
fig_num = 12;
figure(fig_num)
clf;
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')
    
horig(1) = fcn_Points_plotSetsXY(origXYdatasets(1),figure(fig_num))
set(horig(1),'color','blue');
hcorr(1) = fcn_Points_plotSetsXY(biasedXYdataset,figure(fig_num))
set(hcorr(1),'marker','*')
set(hcorr(1),'color','blue');
legend([horig(1) hcorr(1)],{'Original data','Biased data'})
title('fcn_Points_adjustPointSetStatistics: Plot of the original and biased datasets','Interpreter','none');

%% This is a script to exercise the function, fcn_Points_addRadialNoise 

% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;

% Example 1
% shows addition of radial noise

fig_num = 14;
new_set1 = {[0.0 0.0]};


RadiusMaxNoise = 0.5;
for i =  1:100
    datasets_out = fcn_Points_addRadialNoise(new_set1,RadiusMaxNoise);
    fcn_Points_plotSetsXY(new_set1,fig_num);
    fcn_Points_plotSetsXY(datasets_out,fig_num);
    title('fcn_Points_addRadialNoise: Demonstrates the addition of radial noise','Interpreter','none');
    axis([-1 1 -1 1])
    axis square
end

% Example 2
% Addition of radial noise for one data set

fig_num = 15;
RadiusMaxNoise = 1;
new_set2 = {origXYdatasets{1}}; %#ok<CCAT1> 
datasets_out2 = fcn_Points_addRadialNoise(new_set2,RadiusMaxNoise);
fcn_Points_plotSetsXY(new_set2,fig_num); 
fcn_Points_plotSetsXY(datasets_out2,fig_num);
title('fcn_Points_addRadialNoise: Adding noise to the one of the datasets','Interpreter','none')

% Example 3
% Addition of radial noise for one data set
fig_num = 16;
RadiusMaxNoise = 1;

new_set3 = origXYdatasets;
datasets_out3 = fcn_Points_addRadialNoise(new_set3,RadiusMaxNoise);
fcn_Points_plotSetsXY(new_set3,fig_num);
fcn_Points_plotSetsXY(datasets_out3,fig_num);
title('fcn_Points_addRadialNoise: Adding noise to all the datasets','Interpreter','none')


%% This is a script to exercise the function, fcn_Points_calcPairStatistics

% Static test data for testing
% CalcpairedXYdata = [5.93317972350230,88.8321167883212,6.51881720430108,88.6921965317919;10.5414746543779,91.7518248175182,10.3718637992832,90.8598265895954;8.00691244239631,85.0364963503650,8.40053763440860,85.0794797687861;11.9239631336406,86.7883211678832,11.6263440860215,86.2355491329480;9.61981566820277,80.6569343065693,10.1030465949821,81.0332369942196;14.4585253456221,82.1167883211679,13.9560931899642,81.6112716763006;11.6935483870968,75.4014598540146,12.2535842293907,75.5419075144509;16.5322580645161,76.8613138686131,16.0170250896057,76.9869942196532;14.2281105990783,70.4379562043796,14.6729390681004,71.2066473988439;18.8364055299539,73.0656934306569,18.4363799283154,72.5072254335260;16.7626728110599,65.4744525547445,17.1818996415771,66.2933526011560;21.6013824884793,68.1021897810219,21.3037634408602,67.1604046242775;19.0668202764977,61.3868613138686,19.6908602150538,61.9580924855491;21.6013824884793,57.2992700729927,21.9310035842294,57.9118497109826;26.9009216589862,59.0510948905109,26.3216845878136,58.4898843930636;29.4354838709677,56.1313868613139,29.0098566308244,55.5997109826590;28.0529953917051,50.0000000000000,28.6514336917563,49.8193641618497;32.4308755760369,53.5036496350365,31.6980286738351,53.1430635838150;36.1175115207373,51.7518248175182,35.8198924731183,50.9754335260115;34.9654377880184,46.7883211678832,34.5654121863799,47.3627167630058;52.9377880184332,50.2919708029197,52.3073476702509,50.8309248554913;53.1682027649770,57.8832116788321,52.9345878136201,57.0447976878613;55.7027649769585,52.6277372262774,55.1747311827957,52.5650289017341;55.7027649769585,60.8029197080292,55.6227598566308,59.9349710982659;58.9285714285714,55.8394160583942,58.2213261648746,56.4667630057803;58.9285714285714,63.4306569343066,58.6693548387097,62.6806358381503;62.1543778801843,57.8832116788321,61.5367383512545,57.7673410404624;65.3801843317972,58.4671532846715,64.6729390681004,58.6343930635838;72.0622119815668,55.2554744525547,71.8413978494624,55.7442196531792;77.3617511520737,60.5109489051095,76.6801075268817,59.9349710982659;80.3571428571429,56.7153284671533,79.9059139784946,56.0332369942196;83.3525345622120,48.2481751824817,82.6836917562724,48.5187861271676;88.4216589861751,52.6277372262774,87.8808243727599,51.8424855491329;87.0391705069124,45.9124087591241,86.6263440860215,46.0621387283237;91.6474654377880,50.0000000000000,91.1066308243728,49.5303468208092;91.4170506912442,44.7445255474453,90.8378136200717,44.6170520231214;95.5645161290323,48.2481751824817,95.0492831541219,47.6517341040462;95.1036866359447,42.9927007299270,94.6012544802867,43.3164739884393];

% Plotting a sine wave
t=0:0.01:1; 
f=1; 
x=sin(2*pi*f*t);

figure(88); 
plot(t,x,'LineWidth',1.5);
title('Sinusoidal wave with a frequency of 1 Hz','Interpreter','none');

% creating a dataset using t and x values 
data_set = [t;x]';

% chainging the data set into cell array
XYdata = {data_set};

% Adding bias to XYdata
biasx = -2.3;
biasy = 0.9;

noiseMean = zeros(1,2);
noiseVariance = zeros(1,2);

biasedXYdata = fcn_Points_adjustPointSetStatistics(XYdata,[biasx biasy],noiseMean,noiseVariance);


% Plotting the original and biased data on the same figure
fig_num = 21;
figure(fig_num)
axis equal
grid on
xlabel('t')
ylabel('x')

Orig = fcn_Points_plotSetsXY(XYdata,figure(fig_num));
set(Orig,'color','blue');
biaS = fcn_Points_plotSetsXY(biasedXYdata,figure(fig_num));
set(biaS,'marker','*')
set(biaS,'color','blue');

legend([Orig biaS],{'Original data','Biased data'})
title('fcn_Points_adjustPointSetStatistics: Plot of the original and biased datasets','Interpreter','none');


% Run the statistics on the original and biased data

OrignBiasData = [XYdata{1},biasedXYdata{1}];

[errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(OrignBiasData);

fprintf(1,'\n The RMS error is %.4f. \n',errRMS);
fprintf(1,'The variance is %.4f. \n',errVar);
fprintf(1,'The meanShift is %.4f. \n',meanShift);

% Noisy data
RadiusMaxNoise = 0.5;
NoisyXYData = fcn_Points_addRadialNoise(XYdata,RadiusMaxNoise);

axis equal
grid on
fig_num = 51;    
Orig = fcn_Points_plotSetsXY(XYdata,fig_num);
set(Orig,'color','blue');
noisY = fcn_Points_plotSetsXY(NoisyXYData,fig_num);
set(noisY,'marker','*')
set(noisY,'color','red')

legend([Orig noisY],{'Original data','Noisy data'})
title('fcn_Points_addRadialNoise: Plot of the original and noisy datasets','Interpreter','none');

OrignNoisyData = [XYdata{1},NoisyXYData{1}];

[errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(OrignNoisyData);

fprintf(1,'\n The RMS error is %.4f. \n',errRMS);
fprintf(1,'The variance is %.4f. \n',errVar);
fprintf(1,'The meanShift is %.4f. \n',meanShift);



% % Post-process the shift to extract a distance and an angle for easy visual
% % checking of the data
% shiftDist = norm(meanShift,2);
% shiftAngle = atan2(meanShift(2),meanShift(1));

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
figure(11);
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

title('fcn_Points_pairXYdata: Pairing the nearest points','Interpreter','none')

%% This is a script to exercise the function, fcn_Points_pairXYdata for different pairRadius

% Load up some data (simple xy points for now)
load testDatasetVehicle.mat

% Call the pairing function to obtain the matrix of paired XY data without
% a limiting radius
%[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyData{1},xyData{2});

% Define a maximum radius within which to consider data points as pairs
pairRadius = 3.0;

% Call the pairing function to obtain the matrix of paired XY data
[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyData{1},xyData{2},pairRadius);


% Now plot the matches with stars in various colors (to identify the pairs)
figure(11);
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

title('fcn_Points_pairXYdata: Pairing the nearest points for larger pairRadius','Interpreter','none')

%% This is a script to exercise the function, fcn_Points_calcPairStatistics for paired data 

% Static test data for testing
%  CalcpairedXYdata = [5.93317972350230,88.8321167883212,6.51881720430108,88.6921965317919;10.5414746543779,91.7518248175182,10.3718637992832,90.8598265895954;8.00691244239631,85.0364963503650,8.40053763440860,85.0794797687861;11.9239631336406,86.7883211678832,11.6263440860215,86.2355491329480;9.61981566820277,80.6569343065693,10.1030465949821,81.0332369942196;14.4585253456221,82.1167883211679,13.9560931899642,81.6112716763006;11.6935483870968,75.4014598540146,12.2535842293907,75.5419075144509;16.5322580645161,76.8613138686131,16.0170250896057,76.9869942196532;14.2281105990783,70.4379562043796,14.6729390681004,71.2066473988439;18.8364055299539,73.0656934306569,18.4363799283154,72.5072254335260;16.7626728110599,65.4744525547445,17.1818996415771,66.2933526011560;21.6013824884793,68.1021897810219,21.3037634408602,67.1604046242775;19.0668202764977,61.3868613138686,19.6908602150538,61.9580924855491;21.6013824884793,57.2992700729927,21.9310035842294,57.9118497109826;26.9009216589862,59.0510948905109,26.3216845878136,58.4898843930636;29.4354838709677,56.1313868613139,29.0098566308244,55.5997109826590;28.0529953917051,50.0000000000000,28.6514336917563,49.8193641618497;32.4308755760369,53.5036496350365,31.6980286738351,53.1430635838150;36.1175115207373,51.7518248175182,35.8198924731183,50.9754335260115;34.9654377880184,46.7883211678832,34.5654121863799,47.3627167630058;52.9377880184332,50.2919708029197,52.3073476702509,50.8309248554913;53.1682027649770,57.8832116788321,52.9345878136201,57.0447976878613;55.7027649769585,52.6277372262774,55.1747311827957,52.5650289017341;55.7027649769585,60.8029197080292,55.6227598566308,59.9349710982659;58.9285714285714,55.8394160583942,58.2213261648746,56.4667630057803;58.9285714285714,63.4306569343066,58.6693548387097,62.6806358381503;62.1543778801843,57.8832116788321,61.5367383512545,57.7673410404624;65.3801843317972,58.4671532846715,64.6729390681004,58.6343930635838;72.0622119815668,55.2554744525547,71.8413978494624,55.7442196531792;77.3617511520737,60.5109489051095,76.6801075268817,59.9349710982659;80.3571428571429,56.7153284671533,79.9059139784946,56.0332369942196;83.3525345622120,48.2481751824817,82.6836917562724,48.5187861271676;88.4216589861751,52.6277372262774,87.8808243727599,51.8424855491329;87.0391705069124,45.9124087591241,86.6263440860215,46.0621387283237;91.6474654377880,50.0000000000000,91.1066308243728,49.5303468208092;91.4170506912442,44.7445255474453,90.8378136200717,44.6170520231214;95.5645161290323,48.2481751824817,95.0492831541219,47.6517341040462;95.1036866359447,42.9927007299270,94.6012544802867,43.3164739884393];

CalcpairedXYdata = pairedXYdata(1:numMatches,:);

% Run the statistics on the paired data
[errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(CalcpairedXYdata);

fprintf(1,'\n The RMS error is %.4f. \n',errRMS);
fprintf(1,'The variance is %.4f. \n',errVar);
fprintf(1,'The meanShift is %.4f. \n',meanShift);


% Post-process the shift to extract a distance and an angle for easy visual
% checking of the data
shiftDist = norm(meanShift,2);
shiftAngle = atan2(meanShift(2),meanShift(1));

%% This is a script to exercise the function, fcn_Points_pairXYdata with Bias - Make Subplots

% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;
biasx = -1.2;
biasy = 0.6;

noiseMean = zeros(1,2);
noiseVariance = zeros(1,2);
%Corrupt one set of test data with only a systematic bias (for each)
biasedXYdataset = fcn_Points_adjustPointSetStatistics(origXYdatasets(1),[biasx biasy],noiseMean,noiseVariance);


% Call the pairing function to obtain pair the original data with the
% biased data
[pairedXYdataBias, numMatchesBias, nonMatchesABias, nonMatchesBBias] = fcn_Points_pairXYdata(origXYdatasets{1},biasedXYdataset{:});

% Calculate the statistics for the biased data set relative to the original
[errRMSBias,errVarBias,meanShiftBias] = fcn_Points_calcPairStatistics(pairedXYdataBias(1:numMatchesBias,:));



% fh(1) = figure;
fig_num = 12;
figure(fig_num)
set(gcf,'Position',[100 100 1200 500])
clf;
% Plot to provide a visual inspection of the bias
subplot(1,2,1)

grid on
xlabel('x [m]')
ylabel('y [m]')
    
horig(1) = fcn_Points_plotSetsXY(origXYdatasets(1),figure(fig_num));
set(horig(1),'color','blue');
hcorr(1) = fcn_Points_plotSetsXY(biasedXYdataset,figure(fig_num));
set(hcorr(1),'marker','*')
set(hcorr(1),'color','blue');
legend([horig(1) hcorr(1)],{'Original data','Biased data'})

subplot(1,2,2)
% Plotting original data paired with biased data
% figure(13);
% clf

grid on
hold on
xlabel('x [m]')
ylabel('y [m]')
quiver(pairedXYdataBias(:,1),pairedXYdataBias(:,2),pairedXYdataBias(:,3)-pairedXYdataBias(:,1),pairedXYdataBias(:,4)-pairedXYdataBias(:,2),0,'k','linewidth',1)
for i = 1:size(pairedXYdataBias,1)
    plot([pairedXYdataBias(i,1) pairedXYdataBias(i,3)],[pairedXYdataBias(i,2) pairedXYdataBias(i,4)],'*')
end

sgtitle('fcn_Points_pairXYdata with Bias: Plotting the biased dataset and pairing the biased datatset with original','Interpreter','none')

%% This script simulates a vehicle following a trajectory by detecting the nearest points usings KNN algoorithm

% Available in the functions folder
% Change the ROI (Region of interest) around the rover. 

script_AssociatePointsIncremental 


%%%%%%%%%%%%%%%%%%%%%% END OF SCRIPT %%%%%%%%%%%%%%%%%%%%%%%







%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
% Clear out the variables
clear global flag* FLAG*
clear flag*
clear path

% Clear out any path directories under Utilities
path_dirs = regexp(path,'[;]','split');
utilities_dir = fullfile(pwd,filesep,'Utilities');
for ith_dir = 1:length(path_dirs)
    utility_flag = strfind(path_dirs{ith_dir},utilities_dir);
    if ~isempty(utility_flag)
        rmpath(path_dirs{ith_dir});
    end
end

% Delete the Utilities folder, to be extra clean!
if  exist(utilities_dir,'dir')
    [status,message,message_ID] = rmdir(utilities_dir,'s');
    if 0==status
        error('Unable remove directory: %s \nReason message: %s \nand message_ID: %s\n',utilities_dir, message,message_ID);
    end
end

end % Ends fcn_INTERNAL_clearUtilitiesFromPathAndFolders

%% fcn_INTERNAL_initializeUtilities
function  fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders)
% Reset all flags for installs to empty
clear global FLAG*

fprintf(1,'Installing utilities necessary for code ...\n');

% Dependencies and Setup of the Code
% This code depends on several other libraries of codes that contain
% commonly used functions. We check to see if these libraries are installed
% into our "Utilities" folder, and if not, we install them and then set a
% flag to not install them again.

% Set up libraries
for ith_library = 1:length(library_name)
    dependency_name = library_name{ith_library};
    dependency_subfolders = library_folders{ith_library};
    dependency_url = library_url{ith_library};

    fprintf(1,'\tAdding library: %s ...',dependency_name);
    fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
    clear dependency_name dependency_subfolders dependency_url
    fprintf(1,'Done.\n');
end

% Set dependencies for this project specifically
fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders);

disp('Done setting up libraries, adding each to MATLAB path, and adding current repo folders to path.');
end % Ends fcn_INTERNAL_initializeUtilities


function fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url, varargin)
%% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES - MATLAB package installer from URL
%
% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES installs code packages that are
% specified by a URL pointing to a zip file into a default local subfolder,
% "Utilities", under the root folder. It also adds either the package
% subfoder or any specified sub-subfolders to the MATLAB path.
%
% If the Utilities folder does not exist, it is created.
% 
% If the specified code package folder and all subfolders already exist,
% the package is not installed. Otherwise, the folders are created as
% needed, and the package is installed.
% 
% If one does not wish to put these codes in different directories, the
% function can be easily modified with strings specifying the
% desired install location.
% 
% For path creation, if the "DebugTools" package is being installed, the
% code installs the package, then shifts temporarily into the package to
% complete the path definitions for MATLAB. If the DebugTools is not
% already installed, an error is thrown as these tools are needed for the
% path creation.
% 
% Finally, the code sets a global flag to indicate that the folders are
% initialized so that, in this session, if the code is called again the
% folders will not be installed. This global flag can be overwritten by an
% optional flag input.
%
% FORMAT:
%
%      fcn_DebugTools_installDependencies(...
%           dependency_name, ...
%           dependency_subfolders, ...
%           dependency_url)
%
% INPUTS:
%
%      dependency_name: the name given to the subfolder in the Utilities
%      directory for the package install
%
%      dependency_subfolders: in addition to the package subfoder, a list
%      of any specified sub-subfolders to the MATLAB path. Leave blank to
%      add only the package subfolder to the path. See the example below.
%
%      dependency_url: the URL pointing to the code package.
%
%      (OPTIONAL INPUTS)
%      flag_force_creation: if any value other than zero, forces the
%      install to occur even if the global flag is set.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      This code will automatically get dependent files from the internet,
%      but of course this requires an internet connection. If the
%      DebugTools are being installed, it does not require any other
%      functions. But for other packages, it uses the following from the
%      DebugTools library: fcn_DebugTools_addSubdirectoriesToPath
%
% EXAMPLES:
%
% % Define the name of subfolder to be created in "Utilities" subfolder
% dependency_name = 'DebugTools_v2023_01_18';
%
% % Define sub-subfolders that are in the code package that also need to be
% % added to the MATLAB path after install; the package install subfolder
% % is NOT added to path. OR: Leave empty ({}) to only add 
% % the subfolder path without any sub-subfolder path additions. 
% dependency_subfolders = {'Functions','Data'};
%
% % Define a universal resource locator (URL) pointing to the zip file to
% % install. For example, here is the zip file location to the Debugtools
% % package on GitHub:
% dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_18.zip?raw=true';
%
% % Call the function to do the install
% fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
%
% This function was written on 2023_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_23:
% -- wrote the code originally
% 2023_04_20:
% -- improved error handling
% -- fixes nested installs automatically

% TO DO
% -- Add input argument checking

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,4);
end

%% Set the global variable - need this for input checking
% Create a variable name for our flag. Stylistically, global variables are
% usually all caps.
flag_varname = upper(cat(2,'flag_',dependency_name,'_Folders_Initialized'));

% Make the variable global
eval(sprintf('global %s',flag_varname));

if nargin==4
    if varargin{1}
        eval(sprintf('clear global %s',flag_varname));
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist(flag_varname,'var') || isempty(eval(flag_varname))
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;

    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');

        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        end

    end

    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    if ~exist(dependency_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);

        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        end

    end

    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders{1})
        flag_allFoldersThere = 0;
    else
        for ith_folder = 1:length(dependency_subfolders)
            subfolder_name = dependency_subfolders{ith_folder};
            
            % Create the entire path
            subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
            
            % Check if the folder and file exists that is typically created when
            % unzipping.
            if ~exist(subfunction_folder,'dir')
                flag_allFoldersThere = 0;
            end
        end
    end

    % Do we need to unzip the files?
    if flag_allFoldersThere==0
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
        % CANT GET THIS TO WORK --> unzip(zip_file_url, debugTools_folder_name);

        % Is the file there?
        if ~exist(zip_file_name,'file')
            error(['The zip file: %s for dependency: %s did not download correctly.\n' ...
                'This is usually because permissions are restricted on ' ...
                'the current directory. Check the code install ' ...
                '(see README.md) and try again.\n'],zip_file_name, dependency_name);
        end

        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);

        % Did this work? If so, directory should not be empty
        directory_contents = dir(dependency_folder_name);
        if isempty(directory_contents)
            error(['The necessary dependency: %s has an error in install ' ...
                'where the zip file downloaded correctly, ' ...
                'but the unzip operation did not put any content ' ...
                'into the correct folder. ' ...
                'This suggests a bad zip file or permissions error ' ...
                'on the local computer.\n'],dependency_name);
        end

        % Check if is a nested install (for example, installing a folder
        % "Toolsets" under a folder called "Toolsets"). This can be found
        % if there's a folder whose name contains the dependency_name
        flag_is_nested_install = 0;
        for ith_entry = 1:length(directory_contents)
            if contains(directory_contents(ith_entry).name,dependency_name)
                if directory_contents(ith_entry).isdir
                    flag_is_nested_install = 1;
                    install_directory_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name);
                    install_files_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name,'*'); % BUG FIX - For Macs, must be *, not *.*
                    install_location_to = fullfile(directory_contents(ith_entry).folder);
                end
            end
        end

        if flag_is_nested_install
            [status,message,message_ID] = movefile(install_files_from,install_location_to);
            if 0==status
                error(['Unable to move files from directory: %s\n ' ...
                    'To: %s \n' ...
                    'Reason message: %s\n' ...
                    'And message_ID: %s\n'],install_files_from,install_location_to, message,message_ID);
            end
            [status,message,message_ID] = rmdir(install_directory_from);
            if 0==status
                error(['Unable remove directory: %s \n' ...
                    'Reason message: %s \n' ...
                    'And message_ID: %s\n'],install_directory_from,message,message_ID);
            end
        end

        % Make sure the subfolders were created
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders{1})
            for ith_folder = 1:length(dependency_subfolders)
                subfolder_name = dependency_subfolders{ith_folder};
                
                % Create the entire path
                subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
                
                % Check if the folder and file exists that is typically created when
                % unzipping.
                if ~exist(subfunction_folder,'dir')
                    flag_allFoldersThere = 0;
                end
            end
        end
         % If any are not there, then throw an error
        if flag_allFoldersThere==0
            error(['The necessary dependency: %s has an error in install, ' ...
                'or error performing an unzip operation. The subfolders ' ...
                'requested by the code were not found after the unzip ' ...
                'operation. This suggests a bad zip file, or a permissions ' ...
                'error on the local computer, or that folders are ' ...
                'specified that are not present on the remote code ' ...
                'repository.\n'],dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
        end

    end


    % For path creation, if the "DebugTools" package is being installed, the
    % code installs the package, then shifts temporarily into the package to
    % complete the path definitions for MATLAB. If the DebugTools is not
    % already installed, an error is thrown as these tools are needed for the
    % path creation.
    %
    % In other words: DebugTools is a special case because folders not
    % added yet, and we use DebugTools for adding the other directories
    if strcmp(dependency_name(1:10),'DebugTools')
        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');

        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error(['Package installer requires DebugTools package to be ' ...
                'installed first. Please install that before ' ...
                'installing this package']);
        end
    end


    % Finally, the code sets a global flag to indicate that the folders are
    % initialized.  Check this using a command "exist", which takes a
    % character string (the name inside the '' marks, and a type string -
    % in this case 'var') and checks if a variable ('var') exists in matlab
    % that has the same name as the string. The ~ in front of exist says to
    % do the opposite. So the following command basically means: if the
    % variable named 'flag_CodeX_Folders_Initialized' does NOT exist in the
    % workspace, run the code in the if statement. If we look at the bottom
    % of the if statement, we fill in that variable. That way, the next
    % time the code is run - assuming the if statement ran to the end -
    % this section of code will NOT be run twice.

    eval(sprintf('%s = 1;',flag_varname));
end


%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots

    % Nothing to do!



end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function fcn_DebugTools_installDependencies

