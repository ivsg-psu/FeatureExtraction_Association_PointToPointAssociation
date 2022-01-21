% script_test_fcn_Dataset_pairXYdata.m
% This is a script to exercise the function:
% fcn_Dataset_pairXYdata.m

% This script was written on 2022_01_21 by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history
%     2022_01_21
%     -- first write of script

% Clear out the workspace
clearvars

% Load up some data (simple xy points for now)
% %xyData = fcn_Dataset_fillSampleSets;
%load testDataset1.mat
load testDatasetVehicle.mat

% Call the pairing function to obtain the matrix of paired XY data without
% a limiting radius
%[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Dataset_pairXYdata(xyData{1},xyData{2});

% Define a maximum radius within which to consider data points as pairs
pairRadius = 1.0;

% Call the pairing function to obtain the matrix of paired XY data
[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Dataset_pairXYdata(xyData{1},xyData{2},pairRadius);

%% Run some statistics on the paired data
[errRMS,errVar,meanShift] = fcn_Dataset_calcPairStatistics(pairedXYdata(1:numMatches,:));

shiftDist = norm(meanShift,2);
shiftAngle = atan2(meanShift(2),meanShift(1));

%% Now plot the matches with stars in various colors (to identify the pairs)
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