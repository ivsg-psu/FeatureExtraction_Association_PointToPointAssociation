% script_test_fcn_Points_adjustPointSetStatistics.m
% This is a script to exercise the function:
% fcn_Points_adjustPointSetStatistics.m

% This script was written on 2022_01_28 by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history
%     2022_01_28
%     -- wrote the code

% Clear out the workspace
close all
clearvars

% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;

%% Corrupt one set of test data with only a systematic bias (for each)
biasedXYdataset = fcn_Points_adjustPointSetStatistics(origXYdatasets(1),[-0.1 0.1],zeros(1,2),zeros(1,2));

% fh(1) = figure;
% axis equal
% grid on
% xlabel('x [m]')
% ylabel('y [m]')
%     
% horig(1) = fcn_Points_plotSetsXY(origXYdatasets(1),fh(1))
% set(horig(1),'color','blue');
% hcorr(1) = fcn_Points_plotSetsXY(biasedXYdataset,fh(1))
% set(hcorr(1),'marker','*')
% set(hcorr(1),'color','blue');
% legend([horig(1) hcorr(1)],{'Original data','Biased data'})

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

%% Corrupt one set of test data with Gaussian noise in one direction (x-axis here)
noise1DXYdataset = fcn_Points_adjustPointSetStatistics(origXYdatasets(2),zeros(1,2),zeros(1,2),[0.5 0]);

% fh(2) = figure;
% axis equal
% grid on
% xlabel('x [m]')
% ylabel('y [m]')
%     
% horig(2) = fcn_Points_plotSetsXY(origXYdatasets(2),fh(2))
% set(horig(2),'color','blue');
% hcorr(2) = fcn_Points_plotSetsXY(noise1DXYdataset,fh(2))
% set(hcorr(2),'marker','*')
% set(hcorr(2),'color','blue');
% legend([horig(2) hcorr(2)],{'Original data','Data with 1D (x) noise'})

% Call the pairing function to obtain pair the original data with the
% biased data
[pairedXYdataNoise1D, numMatchesNoise1D, nonMatchesANoise1D, nonMatchesBNoise1D] = fcn_Points_pairXYdata(origXYdatasets{2},noise1DXYdataset{:});

% Calculate the statistics for the biased data set relative to the original
[errRMSNoise1D,errVarNoise1D,meanShiftNoise1D] = fcn_Points_calcPairStatistics(pairedXYdataNoise1D(1:numMatchesNoise1D,:));

% Plot to provide a visual inspection of the 1D noise
figure(18);
clf
hold on
quiver(pairedXYdataNoise1D(:,1),pairedXYdataNoise1D(:,2),pairedXYdataNoise1D(:,3)-pairedXYdataNoise1D(:,1),pairedXYdataNoise1D(:,4)-pairedXYdataNoise1D(:,2),0,'k','linewidth',1)
for i = 1:size(pairedXYdataNoise1D,1)
    plot([pairedXYdataNoise1D(i,1) pairedXYdataNoise1D(i,3)],[pairedXYdataNoise1D(i,2) pairedXYdataNoise1D(i,4)],'*')
end

%% Corrupt one set of test data with 2-dimensional Gaussian noise
noise2DXYdataset = fcn_Points_adjustPointSetStatistics(origXYdatasets(2),zeros(1,2),zeros(1,2),[0.1 0.2]);

% fh(3) = figure;
% axis equal
% grid on
% xlabel('x [m]')
% ylabel('y [m]')
%     
% horig(3) = fcn_Points_plotSetsXY(origXYdatasets(2),fh(3))
% set(horig(3),'color','blue');
% hcorr(3) = fcn_Points_plotSetsXY(noise2DXYdataset,fh(3))
% set(hcorr(3),'marker','*')
% set(hcorr(3),'color','blue');
% legend([horig(3) hcorr(3)],{'Original data','Data with 2D (x) noise'})

% Call the pairing function to obtain pair the original data with the
% biased data
[pairedXYdataNoise2D, numMatchesNoise2D, nonMatchesANoise2D, nonMatchesBNoise2D] = fcn_Points_pairXYdata(origXYdatasets{2},noise2DXYdataset{:});

% Calculate the statistics for the biased data set relative to the original
[errRMSNoise2D,errVarNoise2D,meanShiftNoise2D] = fcn_Points_calcPairStatistics(pairedXYdataNoise2D(1:numMatchesNoise2D,:));

% Plot to provide a visual inspection of the 2D noise
figure(19);
clf
hold on
quiver(pairedXYdataNoise2D(:,1),pairedXYdataNoise2D(:,2),pairedXYdataNoise2D(:,3)-pairedXYdataNoise2D(:,1),pairedXYdataNoise2D(:,4)-pairedXYdataNoise2D(:,2),0,'k','linewidth',1)
for i = 1:size(pairedXYdataNoise2D,1)
    plot([pairedXYdataNoise2D(i,1) pairedXYdataNoise2D(i,3)],[pairedXYdataNoise2D(i,2) pairedXYdataNoise2D(i,4)],'*')
end